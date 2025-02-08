// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testModels.C

  \brief Hard-coded basic FE models for various unit testing.

  \details This file contains a set of classes, each one representing a simple
  FE model which is suitable for unit testing of various FE model handling,
  without the need of reading an input file. The FE models are created by the
  respective constructors of the FFlLinkHandler sub-classes.

  \author Knut Morten Okstad

  \date 16 Aug 2021
*/

#include "FFaLib/FFaOS/FFaFortran.H"
#include "FFaLib/FFaAlgebra/FFaVec3.H"
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"

#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlShells.H"
#include "FFlLib/FFlFEParts/FFlRGD.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"
#include "FFlLib/FFlFEParts/FFlPTHICK.H"

#ifdef FF_NAMESPACE
using namespace FF_NAMESPACE;
#endif


/*!
  \brief Static helper generating a quadrilateral element in a regular mesh.
*/

static FFlElementBase* newQuadElement (int e, int i, int j, int n1,
                                       bool parabolic = false)
{
  FFlElementBase* elm = NULL;
  if (parabolic)
  {
    elm = new FFlQUAD8(e);
    elm->setNode(1,n1*(i+0)+j+1);
    elm->setNode(2,n1*(i+0)+j+2);
    elm->setNode(3,n1*(i+0)+j+3);
    elm->setNode(4,n1*(i+1)+j+3);
    elm->setNode(5,n1*(i+2)+j+3);
    elm->setNode(6,n1*(i+2)+j+2);
    elm->setNode(7,n1*(i+2)+j+1);
    elm->setNode(8,n1*(i+1)+j+1);
  }
  else
  {
    elm = new FFlQUAD4(e);
    elm->setNode(1,n1*(i+0)+j+1);
    elm->setNode(2,n1*(i+0)+j+2);
    elm->setNode(3,n1*(i+1)+j+2);
    elm->setNode(4,n1*(i+1)+j+1);
  }

  return elm;
}


/*!
  \brief Creates an FE model of a shell strip with rigid ends.
*/

class ShellStrip : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] L Total length of the strip
  //! \param[in] b Width of the strip
  //! \param[in] t Shell thickness
  //! \param[in] nL Number of elements in length direction
  //! \param[in] nb Number of elements in width direction
  ShellStrip(double L, double b, double t, int nL, int nb)
  {
    std::cout <<"   * Generating a plate: L="<< L <<" b="<< b <<" t="<< t
              <<" nL="<< nL <<" nb="<< nb << std::endl;

    FFlPMAT* pmat = new FFlPMAT(1);
    pmat->youngsModule = 2.1e11;
    pmat->poissonsRatio = 0.3;
    pmat->materialDensity = 7850.0;
    this->addAttribute(pmat);

    FFlPTHICK* pthk = new FFlPTHICK(1);
    pthk->thickness = t;
    this->addAttribute(pthk);

    int nnb = nb+1;
    int nnL = nL+1;

    int i, j, e;
    FFlNode* nn = NULL;
    for (j = 0; j <= nb; j++)
      for (i = 0; i <= nL; i++)
      {
        FaVec3 X((double)i*L/(double)nL,(double)j*b/(double)nb,0.0);
        this->addNode(nn = new FFlNode(nnL*j+i+1,X));
      }

    FFlElementBase* elm = NULL;
    for (i = e = 0; i < nb; i++)
      for (j = 0; j < nL; j++)
      {
        elm = newQuadElement(++e,i,j,nnL);
        elm->setAttribute(pthk);
        elm->setAttribute(pmat);
        this->addElement(elm);
      }

    this->addNode(nn = new FFlNode(nnb*nnL+1,FaVec3(0.0,0.5*b,0.0)));
    elm = new FFlRGD(++e);
    elm->setNode(1,nn);
    for (j = 0; j < nnb; j++)
      elm->setNode(2+j,1+j*nnL);
    this->addElement(elm);

    this->addNode(nn = new FFlNode(nnb*nnL+2,FaVec3(L,0.5*b,0.0)));
    elm = new FFlRGD(++e);
    elm->setNode(1,nn);
    for (j = 0; j < nnb; j++)
      elm->setNode(2+j,(1+j)*nnL);
    this->addElement(elm);

    this->resolve();
  }
};


void ffl_setLink (FFlLinkHandler* link);

//! \brief Creates an FE model as a FFlLinkHandler object.
SUBROUTINE(make_femodel,MAKE_FEMODEL) (const int& iPart,
                                       const int& ne1, const int& ne2,
                                       const double& L, const double& b,
                                       const double& t, int& ierr)
{
  FFaCmdLineArg::init(0,NULL);
  FFaCmdLineArg::mute = true;
  FFaCmdLineArg::instance()->addOption ("useANDESformulation",true,
                                        "Shell formulation option");
  FFlNode::init();
  FFlQUAD4::init();
  FFlQUAD8::init();
  FFlRGD::init();
  FFlPMAT::init();
  FFlPTHICK::init();

  ierr = 0;
  switch (iPart) {
  case 1: ffl_setLink(new ShellStrip(L,b,t,ne1,ne2)); break;
  default:
    std::cerr <<" *** Non-implemented FE part "<< iPart << std::endl;
    ierr = -iPart;
  }
}
