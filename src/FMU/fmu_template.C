/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file fmu_template.cpp
  \brief FMU implementation for FEDEM
  \details This file implements the subset of the functions declared in the file
  fmi2Functions.h from the FMI2 standard, necessary to run a FEDEM model as an
  FMU in a co-simulation context.
*/

#undef UNICODE

#include "fmi2Functions.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

#if defined _MSC_VER
#include <windows.h>
typedef HINSTANCE LibHandle;
#elif defined __GNUC__
#include <dlfcn.h>
typedef void* LibHandle;
#else
#error "Platform not supported, neither _MSC_VER nor __GNUC__ defined"
#endif

//Macro for printing debug info to terminal
#ifdef FMU_DEBUG
#define DEBUG_STDOUT(x) std::cout << x << std::endl;
#else 
#define DEBUG_STDOUT(x)
#endif

typedef void (*DLPROC)();
typedef int (*DLPROC_INIT)(int,char**,const char*,
                           const double*,const int,
                           const double*,const int,
                           const double*,int*);
typedef int (*DLPROC_DONE)(bool);
typedef int (*DLPROC_GETSTATESIZE)();
typedef bool (*DLPROC_SETEXTFUNC)(int,double);
typedef bool (*DLPROC_SOLVENEXT)(int*);
typedef double (*DLPROC_EVALFUNC)(int,const char*,double,int*);
typedef bool (*DLPROC_SAVETRANSFORMATIONSTATE)(double*,const int);
typedef int (*DLPROC_RESTARTFROMSTATE)(const double*,const int,const int);

DLPROC_INIT solverInit;
DLPROC_DONE solverDone;
DLPROC_GETSTATESIZE getStateSize;
DLPROC_GETSTATESIZE getTransformationStateSize;
DLPROC_SETEXTFUNC setExtFunc;
DLPROC_SOLVENEXT solveNext;
DLPROC_EVALFUNC evalFunc;
DLPROC_SAVETRANSFORMATIONSTATE saveTransformationState;
DLPROC_RESTARTFROMSTATE restartFromState;


namespace {

  enum fmuStateCode
  {
    FMUSTART= 1u << 0,
    FMUEND = 1u << 1,
    FMUINSTANTIATED = 1u << 2,
    FMUINITIALIZATION = 1u << 3,
    FMUSTEPCOMPLETE= 1u << 4,
    FMUSTEPINPROGRESS = 1u << 5,
    FMUSTEPFAILED = 1u << 6,
    FMUSTEPCANCELED = 1u << 7,
    FMUTERMINATED = 1u << 8,
    FMUERROR = 1u << 9,
    FMUFATAL = 1u << 10
  };
  
  //State of internal fmu parameters and variables, as well as solver state for restart.
  struct modelState
  {
    std::string modelIdentifier;
    std::string modelGuid;
    fmi2Integer numReals;
    fmi2Integer numInputs;
    fmi2Integer numOutputs;
    fmi2Integer numParams;
    fmi2Integer numTransforms;
    
    fmi2Integer* fedemInputIndices;
    fmi2Integer* fedemOutputIndices;
    fmi2Integer* fedemTransformIndices;
    
    fmi2Integer solverStateSize;
    fmi2Real* solverState;
    fmi2Integer transformationStateSize;
    fmi2Real* transformationState;
    fmi2Real t;
    fmi2Real* reals;
  };
  
  struct componentInstance
  {
    fmuStateCode stateCode;
    modelState state;
    modelState initialState;
    fmi2String instanceName;
    fmi2Boolean logging = false;
    const fmi2CallbackFunctions *functions;
  };

  void readConfig(componentInstance* comp, const std::string& path)
  {
    DEBUG_STDOUT("configPath: " + path);

    // Lambda function allocating an integer array.
    auto&& allocateInts = [comp](fmi2Integer nwi)
    {
      return (fmi2Integer*)comp->functions->allocateMemory(nwi,sizeof(fmi2Integer));
    };

    // TODO(RunarHR): readConfig: Make this more robust
    std::string line;
    std::ifstream confFile(path);
    std::getline(confFile,comp->initialState.modelIdentifier);
    std::getline(confFile,comp->initialState.modelGuid);
    std::getline(confFile,line);
    comp->initialState.numReals = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numInputs = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numOutputs = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numParams = atoi(line.c_str());
    std::getline(confFile,line);
    comp->initialState.numTransforms = atoi(line.c_str());
    comp->initialState.fedemInputIndices = allocateInts(comp->initialState.numInputs);
    comp->initialState.fedemOutputIndices = allocateInts(comp->initialState.numOutputs);
    comp->initialState.fedemTransformIndices = allocateInts(comp->initialState.numTransforms);
    comp->state.fedemInputIndices = allocateInts(comp->initialState.numInputs);
    comp->state.fedemOutputIndices = allocateInts(comp->initialState.numOutputs);
    comp->state.fedemTransformIndices = allocateInts(comp->initialState.numTransforms);
    
    for(int i=0; i< comp->initialState.numInputs; i++)
    {
      std::getline(confFile,line);
      comp->initialState.fedemInputIndices[i] = atoi(line.c_str());
    }
    for(int i=0; i< comp->initialState.numOutputs; i++)
    {
      std::getline(confFile,line);
      comp->initialState.fedemOutputIndices[i] = atoi(line.c_str());
    }
    for(int i=0; i< comp->initialState.numTransforms; i++)
    {
      std::getline(confFile,line);
      comp->initialState.fedemTransformIndices[i] = atoi(line.c_str());
    }
  }
  
  void copyModelState(modelState* destination, modelState* source)
  {
    //Copy arrays
    memcpy( destination->reals, source->reals, sizeof(fmi2Real)*source->numReals );
    memcpy( destination->solverState, source->solverState, sizeof(fmi2Real)*source->solverStateSize );
    memcpy( destination->transformationState, source->transformationState, sizeof(fmi2Real)*source->transformationStateSize );
    
    //Copy variables
    destination->modelIdentifier = source->modelIdentifier;
    destination->modelGuid = source->modelGuid;
    destination->numReals = source->numReals;
    destination->numInputs = source->numInputs;
    destination->numOutputs = source->numOutputs;
    destination->numParams = source->numParams;
    destination->numTransforms = source->numTransforms;
    
    destination->t = source->t;
    destination->solverStateSize = source->solverStateSize;
    destination->transformationStateSize = source->transformationStateSize;
    
    memcpy( destination->fedemInputIndices, source->fedemInputIndices, sizeof(fmi2Integer)*source->numInputs );
    memcpy( destination->fedemOutputIndices, source->fedemOutputIndices, sizeof(fmi2Integer)*source->numOutputs );
    memcpy( destination->fedemTransformIndices, source->fedemTransformIndices, sizeof(fmi2Integer)*source->numTransforms );
  }

} // closing brace for anonymous namespace


#ifdef __cplusplus
extern "C" {
#endif

  fmi2Component fmi2Instantiate(fmi2String  instanceName, fmi2Type fmuType, fmi2String  GUID, fmi2String  fmuResourceLocation, const fmi2CallbackFunctions* functions, fmi2Boolean visible, fmi2Boolean loggingOn)
  {
    DEBUG_STDOUT("Hello, FMU world!");

    const char* solverPath = getenv("FEDEM_SOLVER");
    if (!solverPath)
    {
      std::cerr <<" *** Environment variable FEDEM_SOLVER not defined."<< std::endl;
      return 0;
    }

    DEBUG_STDOUT("solverPath: " + std::string(solverPath));

    componentInstance* comp = (componentInstance*)functions->allocateMemory(1,sizeof(componentInstance));
    comp->logging = loggingOn;
    comp->functions = functions;
    comp->instanceName = instanceName;

    // Lambda function returning the parent path of given pathname.
    auto&& parentPath = [](const std::string& path)
    {
      size_t pos = path.find_last_of("/\\");
      return path.substr(0, pos == std::string::npos ? 0 : pos);
    };

    // Lambda function appending a filename to a path name.
    auto&& appendPath = [](const std::string& path, const std::string& fname)
    {
#if defined _MSC_VER
      return path + "\\" + fname;
#else
      return path + "/" + fname;
#endif
    };

    //Get resource folder path
#if defined _MSC_VER
    char path[1024];
    HMODULE hm = NULL;
    if (!GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCSTR)&fmi2GetTypesPlatform, &hm))
      if (comp->logging)
        functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
                          "fmi2Instantiate: Could not retrieve path of running module from GetModuleHandleEx.");
    GetModuleFileNameA(hm, path, sizeof(path));
#elif defined __GNUC__
    Dl_info info;
    if (!dladdr((void*)fmi2Instantiate, &info))
    {
      std::cerr <<" *** Could not find location of FMU-library."<< std::endl;
      return 0;
    }
    const char* path = info.dli_fname;
#endif
    DEBUG_STDOUT("FMU-library location: " + std::string(path));

    std::string platformPath(parentPath(path));
    DEBUG_STDOUT("platformPath: " + platformPath);
    std::string binariesPath(parentPath(platformPath));
    DEBUG_STDOUT("binariesPath: " + binariesPath);
    std::string rootPath(parentPath(binariesPath));
    std::string resourcePath(appendPath(rootPath,"resources"));
    std::string workingDir(appendPath(resourcePath,"model"));
    DEBUG_STDOUT("workingDir: " + workingDir);

    readConfig(comp, appendPath(resourcePath,"config.txt"));

    if (strcmp(comp->initialState.modelGuid.c_str(),GUID) != 0)
    {
      std::cerr <<" *** GUIDs do not match: "<< comp->initialState.modelGuid <<" "<< GUID << std::endl;
      if(comp->logging)functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
                                         "fmi2Instantiate: GUIDs does not match.");
      functions->freeMemory(comp);
      return 0;
    }

    // NOTE(RunarHR): Initialization of solver should be done in EnterInitializationMode, but must be done here to get solverState size. 
#if defined _MSC_VER
    std::string solverBinariesPath(parentPath(solverPath));
    SetDllDirectory(solverBinariesPath.c_str());
    LibHandle h_solver = LoadLibrary(solverPath);
    SetDllDirectory(NULL);
#elif defined __GNUC__
    LibHandle h_solver = dlopen(solverPath,RTLD_LAZY);
#endif
    if (!h_solver) {
      std::cerr <<" *** Could not load solver library "<< solverPath << std::endl;
      return 0;
    }

    // Lambda function of optaining a function pointer from the solver library.
    auto&& getFuncAddress = [h_solver](const char* fName)
    {
#if defined _MSC_VER
      DLPROC address = (DLPROC)GetProcAddress(h_solver,fName);
#else
      DLPROC address = (DLPROC)dlsym(h_solver,fName);
#endif
      if (!address)
        std::cerr <<" *** Could not get function address for "<< fName << std::endl;
      return address;
    };

    solverInit = (DLPROC_INIT)getFuncAddress("solverInit");
    solverDone = (DLPROC_DONE)getFuncAddress("solverDone");
    getStateSize = (DLPROC_GETSTATESIZE)getFuncAddress("getStateSize");
    getTransformationStateSize = (DLPROC_GETSTATESIZE)getFuncAddress("getTransformationStateSize");
    setExtFunc = (DLPROC_SETEXTFUNC)getFuncAddress("setExtFunc");
    solveNext = (DLPROC_SOLVENEXT)getFuncAddress("solveNext");
    evalFunc = (DLPROC_EVALFUNC)getFuncAddress("evalFunc");
    saveTransformationState = (DLPROC_SAVETRANSFORMATIONSTATE)getFuncAddress("saveTransformationState");
    restartFromState = (DLPROC_RESTARTFROMSTATE)getFuncAddress("restartFromState");

    // Lambda function adding option files to argvStart based on existance.
    std::vector<char*> argvStart;
    auto&& addOptionFile=[workingDir,&argvStart,appendPath](const char* opt,const char* fileName)
    {
      std::ifstream fs(appendPath(workingDir,fileName).c_str());
      if (fs.good())
      {
        argvStart.push_back(const_cast<char*>(opt));
        argvStart.push_back(const_cast<char*>(fileName));
      }
    };

    argvStart.reserve(9);
    argvStart.push_back(const_cast<char*>("fedem_solver"));
    argvStart.push_back(const_cast<char*>("-cwd"));
    argvStart.push_back(const_cast<char*>(workingDir.c_str()));
    addOptionFile("-fco","fedem_solver.fco");
    addOptionFile("-fop","fedem_solver.fop");
    addOptionFile("-fao","fedem_solver.fao");

    DEBUG_STDOUT("Initializing solver");
    int status = solverInit(argvStart.size(),argvStart.data(),NULL,NULL,0,NULL,0,NULL,NULL);
    if (status < 0)
    {
      if(comp->logging) functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error",
                                          "Could not initialize solver.");
      std::cerr <<" *** Solver failed to initialize ("<< status <<")."<< std::endl;
      return 0;
    }

    // Lambda function allocating a real array.
    auto&& allocateReals = [functions](fmi2Integer nwr)
    {
      return (fmi2Real*)functions->allocateMemory(nwr,sizeof(fmi2Real));
    };

    comp->initialState.solverStateSize = getStateSize();
    comp->initialState.transformationStateSize = getTransformationStateSize();

    comp->initialState.reals = allocateReals(comp->initialState.numReals);
    comp->initialState.solverState = allocateReals(comp->initialState.solverStateSize);
    comp->initialState.transformationState = allocateReals(comp->initialState.transformationStateSize);
    comp->initialState.t = 0;

    comp->state.reals = allocateReals(comp->initialState.numReals);
    comp->state.solverState = allocateReals(comp->initialState.solverStateSize);
    comp->state.transformationState = allocateReals(comp->initialState.transformationStateSize);
    
    copyModelState( &(comp->state), &(comp->initialState) );
    comp->stateCode = fmuStateCode::FMUINSTANTIATED;
    
    return comp;
  }
  
  //Stop simulation
  fmi2Status fmi2Terminate(fmi2Component c)
  {
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED))
    {
      //Close solver
      //solverDone(true);
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2DoStep(fmi2Component c, fmi2Real currentCommunicationPoint, fmi2Real communicationStepSize, fmi2Boolean newStep)
  {
    DEBUG_STDOUT("fmi2DoStep");

    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & fmuStateCode::FMUSTEPCOMPLETE)
    {
      comp->state.t = currentCommunicationPoint;

      //SET INPUT FUNCTIONS IN FEDEM.
      int err = 0;
      for (int i = 0; i < comp->state.numInputs && err == 0; i++)
        if (!setExtFunc(comp->state.fedemInputIndices[i], comp->state.reals[i]))
          --err;

      if (err == 0) solveNext(&err);

      if (err != 0)
      {
        if(comp->logging)comp->functions->logger(comp->functions->componentEnvironment, comp->instanceName, fmi2Error, "error",
                                                 "fmi2DoStep: Solver step failed.");
        comp->stateCode = fmuStateCode::FMUSTEPFAILED;
        return fmi2Error;
      }

      //GET OUTPUT FROM FEDEM
      for (int i = 0; i < comp->state.numOutputs && err == 0; i++)
        comp->state.reals[comp->state.numInputs+i] = evalFunc(comp->state.fedemOutputIndices[i], NULL, -1.0, &err);

      if (err != 0)
      {
        comp->stateCode = fmuStateCode::FMUSTEPFAILED;
        return fmi2Error;
      }     

      //GET TRANSFORMATION OUTPUT FROM FEDEM.
      if(comp->state.numTransforms > 0)
      {
        saveTransformationState(comp->state.transformationState, comp->state.transformationStateSize);

        size_t idat = 3; //advance past timestep data
        //Loop over all objects in transformationState. Size per object is 14.
        fmi2Real* mem = comp->state.transformationState;
        for (int i = 0; i < (comp->state.transformationStateSize-3)/14; i++)
        {
          idat++; //Skip typeId
          fmi2Integer baseID = (fmi2Integer)mem[idat++];

          //Check if baseID of transformation object match any of the
          //ones specified for the FMU.
          for (int j = 0; j < comp->state.numTransforms; j++)
            if (comp->state.fedemTransformIndices[j] == baseID)
            {
              fmi2Integer offset = comp->state.numInputs + comp->state.numOutputs + comp->state.numParams + j*12;
              memcpy(comp->state.reals+offset, mem+idat, sizeof(fmi2Real)*12);
              break;
            }

          idat += 12;
        }
      }

      comp->stateCode = fmuStateCode::FMUSTEPCOMPLETE;
      return fmi2OK;
    }
    
    // TODO(RunarHR): fmi2DoStep: If communicationStepSize is not equal to predefined step-size: return fmi2Error and write log-message
    
    comp->stateCode = fmuStateCode::FMUSTEPFAILED;
    return fmi2Error;
  }
  
  
  void fmi2FreeInstance(fmi2Component c)
  {
    //Clean up componentInstance, free all memory and resources
    
    if(c == 0)
      return;
    
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      //Free state
      comp->functions->freeMemory( comp->state.reals );
      comp->functions->freeMemory( comp->state.solverState );
      comp->functions->freeMemory( comp->state.transformationState );
      comp->functions->freeMemory( comp->state.fedemInputIndices );
      comp->functions->freeMemory( comp->state.fedemOutputIndices );
      
      //Free initial state
      comp->functions->freeMemory( comp->initialState.reals );
      comp->functions->freeMemory( comp->initialState.solverState );
      comp->functions->freeMemory( comp->initialState.transformationState );
      comp->functions->freeMemory( comp->initialState.fedemInputIndices );
      comp->functions->freeMemory( comp->initialState.fedemOutputIndices );
      
      comp->functions->freeMemory( comp );
      
      //Close solver
      solverDone(true);
    }
    
    return;
  }
  
  fmi2Status fmi2SetupExperiment(fmi2Component c, fmi2Boolean toleranceDefined, fmi2Real tolerance, fmi2Real startTime, fmi2Boolean stopTimeDefined, fmi2Real stopTime)
  {
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode == fmuStateCode::FMUINSTANTIATED)
    {
      // TODO(RunarHR): fmi2SetupExperiment: setup start time etc. Called before initialize
      //Setup solver parameters
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2EnterInitializationMode(fmi2Component c) {
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode == fmuStateCode::FMUINSTANTIATED)
    {
      comp->stateCode = fmuStateCode::FMUINITIALIZATION;
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2ExitInitializationMode(fmi2Component c) {
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode == fmuStateCode::FMUINITIALIZATION)
    {
      // TODO(RunarHR): fmi2ExitInitializationMode: If values corresponding to input channels/external functions have been set, update solver. 
      
      comp->stateCode = fmuStateCode::FMUSTEPCOMPLETE;
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  
  //WARNING: This fails with Ansys TwinBuilder. The fmi2FMUstate pointer is not NULL on first call.
  fmi2Status fmi2GetFMUstate (fmi2Component c, fmi2FMUstate* FMUstate) {
    
    /*From Doc: fmi2GetFMUstate makes a copy of the internal FMU state and returns a pointer to this copy
    (FMUstate). If on entry *FMUstate == NULL, a new allocation is required. If *FMUstate !=
    NULL, then *FMUstate points to a previously returned FMUstate that has not been modified
    since. In particular, fmi2FreeFMUstate had not been called with this FMUstate as an argument.
    [Function fmi2GetFMUstate typically reuses the memory of this FMUstate in this case and
    returns the same pointer to it, but with the actual FMUstate.]*/
    
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      modelState* state = (modelState*)FMUstate;
      
      if(state == 0)
      {
        state = (modelState*)comp->functions->allocateMemory(1,sizeof(modelState));
        state->reals = (fmi2Real*)comp->functions->allocateMemory(comp->state.numReals,sizeof(fmi2Real));
        state->solverState = (fmi2Real*)comp->functions->allocateMemory(comp->state.solverStateSize,sizeof(fmi2Real));
        state->transformationState = (fmi2Real*)comp->functions->allocateMemory(comp->state.transformationStateSize,sizeof(fmi2Real));
      }
      
      copyModelState(state, &(comp->state));
      
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
    
  }
  fmi2Status fmi2SetFMUstate (fmi2Component c, fmi2FMUstate FMUstate) {
    /*From Doc: fmi2SetFMUstate copies the content of the previously copied FMUstate back and uses it as
    actual new FMU state. The FMUstate copy does still exist.*/
    
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      modelState* state = (modelState*)FMUstate;
      
      if(state == 0)
      {
        comp->stateCode = fmuStateCode::FMUERROR;
        return fmi2Error;
      }
      
      //Copy input state to stored state
      copyModelState(&(comp->state), state);
      
      //Reset solver.
      restartFromState(comp->state.solverState, comp->state.solverStateSize,0); // NOTE(RunarHR): writeToRDB is 0. No results saved.
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate) {
    /* From Doc: fmi2FreeFMUstate frees all memory and other resources allocated with the fmi2GetFMUstate
    call for this FMUstate. The input argument to this function is the FMUstate to be freed. If a null
    pointer is provided, the call is ignored. The function returns a null pointer in argument FMUstate.
    These functions are only supported by the FMU, if the optional capability flag
    <fmiModelDescription> <ModelExchange / CoSimulation canGetAndSetFMUstate in =
    "true"> in the XML file is explicitly set to true (see sections 3.3.1 and 4.3.1).*/
    
    componentInstance* comp = (componentInstance *)c;
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      modelState* state = (modelState*)FMUstate;
      
      if(state != 0)
      {
        comp->functions->freeMemory(state->solverState);
        comp->functions->freeMemory(state->transformationState);
        comp->functions->freeMemory(state->reals);
        comp->functions->freeMemory(state);
      }
      
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[])
  {
    DEBUG_STDOUT("fmi2SetReal");
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE))
    {
      if(comp->state.numReals == 0)
      {
        comp->stateCode = fmuStateCode::FMUERROR;
        return fmi2Error;
      }

      for(size_t i = 0; i < nvr; i++)
      {
        if((fmi2Integer)vr[i] >= ( comp->state.numInputs + comp->state.numOutputs + comp->state.numParams + comp->state.numTransforms))
        {
          comp->stateCode = fmuStateCode::FMUERROR;
          return fmi2Error;
        }
        comp->state.reals[vr[i]] = value[i];
      }
      
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[])
  {
    DEBUG_STDOUT("fmi2GetReal");
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      DEBUG_STDOUT("stateCode-correct");
      if(comp->state.numReals == 0)
      {
        DEBUG_STDOUT("num reals == 0");
        comp->stateCode = fmuStateCode::FMUERROR;
        return fmi2Error;
      }

      for(size_t i = 0; i < nvr; i++)
      {
        if((fmi2Integer)vr[i] >= ( comp->state.numInputs + comp->state.numOutputs + comp->state.numParams + comp->state.numTransforms*12 ))
        {
          DEBUG_STDOUT("Vr > numInputs + numOutputs + numParams + numTransforms*12");
          comp->stateCode = fmuStateCode::FMUERROR;
          return fmi2Error;
        }
        value[i] = comp->state.reals[vr[i]];
      }
      
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  const char* fmi2GetTypesPlatform() {
    DEBUG_STDOUT("fmi2GetTypesPlatform");
    return "default";
  }
  
  fmi2Status fmi2Reset(fmi2Component c)
  {
    DEBUG_STDOUT("fmi2Reset");
    componentInstance* comp = (componentInstance *)c;
    
    if(comp->stateCode & (fmuStateCode::FMUINSTANTIATED | fmuStateCode::FMUINITIALIZATION | fmuStateCode::FMUSTEPCOMPLETE | fmuStateCode::FMUSTEPFAILED
                          | fmuStateCode::FMUSTEPCANCELED | fmuStateCode::FMUTERMINATED | fmuStateCode::FMUERROR))
    {
      //Copy initialState to state
      memcpy( comp->state.reals, comp->initialState.reals, sizeof(fmi2Real)*(comp->initialState.numReals));
      memcpy( comp->state.solverState, comp->initialState.solverState, sizeof(fmi2Real)*(comp->initialState.solverStateSize));
      memcpy( comp->state.transformationState, comp->initialState.transformationState, sizeof(fmi2Real)*(comp->initialState.transformationStateSize));
      
      //Reset variables
      comp->state.t = comp->initialState.t;
      
      //Reset solver.
      //restartFromState(comp->state.solverState, comp->state.solverStateSize,0); // NOTE(RunarHR): writeToRDB is 0. No results saved.
      
      comp->stateCode = fmuStateCode::FMUINSTANTIATED;
      
      return fmi2OK;
    }
    
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  const char* fmi2GetVersion() {
    return fmi2Version;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////
  //        NOT IMPLEMENTED                           //
  //////////////////////////////////////////////////////

  fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetInteger
    DEBUG_STDOUT("fmiGetInteger");
    return fmi2OK;
  }
  
  
  fmi2Status fmi2GetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetBoolean
    DEBUG_STDOUT("fmiGetBoolean");
    return fmi2OK;
  }
  
  
  fmi2Status fmi2GetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetString
    return fmi2OK;
  }
  
  // Derivatives is not supported. 
  //canInterpolateInputs and MaxOutputDerivativeOrder must be set to false and 0 in element "CoSimulation" in the xml-file
  fmi2Status fmi2SetRealInputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], const fmi2Real value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2SetRealInputDerivatives
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t  nvr, const fmi2Integer order[], fmi2Real value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetRealOutputDerivatives
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2CancelStep(fmi2Component c)
  {
    // TODO(RunarHR): [Optional] Implement fmi2CancelStep
    // Only relevant if doStep runs asynchronously.
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2GetStatus(fmi2Component c, const fmi2StatusKind s, fmi2Status* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2GetRealStatus(fmi2Component c, const fmi2StatusKind s, fmi2Real* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetRealStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2GetIntegerStatus(fmi2Component c, const fmi2StatusKind s, fmi2Integer* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetIntegerStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2GetBooleanStatus(fmi2Component c, const fmi2StatusKind s, fmi2Boolean* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetBooleanStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2GetStringStatus(fmi2Component c, const fmi2StatusKind s, fmi2String* value)
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetStringStatus
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2SetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2SetInteger
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2SetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2SetBoolean
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2SetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2String value[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2SetString
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  
  fmi2Status fmi2SetDebugLogging(fmi2Component c, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[]) 
  {
    // TODO(RunarHR): [Optional] Implement fmi2SetDebugLogging
    /*componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;*/
    return fmi2OK;
  }
  
  fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t *size) {
    // TODO(RunarHR): [Optional] Implement fmi2SerializedFMUstateSize
    
    /*From Doc: fmi2SerializedFMUstateSize returns the size of the byte vector, in order that FMUstate can
    be stored in it. With this information, the environment has to allocate an fmi2Byte vector of the
    required length size.*/
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  fmi2Status fmi2SerializeFMUstate (fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[], size_t size) {
    // TODO(RunarHR): [Optional] Implement fmi2SerializeFMUstate
    
    /*From Doc: fmi2SerializeFMUstate serializes the data which is referenced by pointer FMUstate and
    copies this data in to the byte vector serializedState of length size, that must be provided by
    the environment.*/
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }
  fmi2Status fmi2DeSerializeFMUstate (fmi2Component c, const fmi2Byte serializedState[], size_t size,
                                      fmi2FMUstate* FMUstate) {
    // TODO(RunarHR): [Optional] Implement fmi2DeSerializeFMUstate
    
    /*From Doc: fmi2DeSerializeFMUstate deserializes the byte vector serializedState of length size,
    constructs a copy of the FMU state and returns FMUstate, the pointer to this copy. [The
    simulation is restarted at this state, when calling fmi2SetFMUState with FMUstate.]
    These functions are only supported by the FMU, if the optional capability flags
    canGetAndSetFMUstate and canSerializeFMUstate in
    <fmiModelDescription><ModelExchange / CoSimulation> in the XML file are explicitly set
    to true (see sections 3.3.1 and 4.3.1).*/
    componentInstance* comp = (componentInstance *)c;
    comp->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }

  fmi2Status fmi2GetDirectionalDerivative(fmi2Component c,
                                          const fmi2ValueReference vUnknown_ref[], size_t,
                                          const fmi2ValueReference vKnown_ref[], size_t,
                                          const fmi2Real dvKnown[], fmi2Real dvUnknown[])
  {
    // TODO(RunarHR): [Optional] Implement fmi2GetDirectionalDerivative
    static_cast<componentInstance*>(c)->stateCode = fmuStateCode::FMUERROR;
    return fmi2Error;
  }

#ifdef __cplusplus
} // closing brace for extern "C"
#endif
