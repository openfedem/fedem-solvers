/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file binaryDB.c
  \brief Global C-functions for binary file IO.

  \details This file contains the C-implementation of the following
  Fortran-callable global functions for performing binary file IO,
  based on low-level fwrite/fread calls.

  - binarydbinterface::openBinaryDB
  - binarydbinterface::closeBinary_DB
  - binarydbinterface::deleteDB
  - binarydbinterface::setBufDB
  - binarydbinterface::flushBinaryDB
  - binarydbinterface::writeIntDB
  - binarydbinterface::writeFloatDB
  - binarydbinterface::writeDoubleD4
  - binarydbinterface::writeDoubleD8
  - binarydbinterface::writeCharDB
  - binarydbinterface::writeTagDB
  - binarydbinterface::readIntDB
  - binarydbinterface::readFloatDB
  - binarydbinterface::readDoubleD4
  - binarydbinterface::readDoubleD8
  - binarydbinterface::readCharDB
  - binarydbinterface::readTagDB
  - binarydbinterface::setPositionDB
  - binarydbinterface::getPositionDB
  - binarydbinterface::putCharDB
  - binarydbinterface::copyBinaryDB

  No further documentation of the above functions is provided here.
  They are documented in the binaryDBInterface.f90 file.

  This file also contains some static data containers managing the binary files,
  and some private functions serving as utilities for the global functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#if !defined(win32) && !defined(win64)
#include <unistd.h>
#endif
#include "FFaLib/FFaOS/FFaTag_C.h"
#include "FFaLib/FFaOS/FFaFortran.H"

#if _MSC_VER > 1310
#define unlink  _unlink
#define tempnam _tempnam
#endif

#ifndef _MAX_DBFIL
/*! Max number of binary files that can be open at the same time */
#define _MAX_DBFIL 50
#endif

static char* fname[_MAX_DBFIL]; /*!< Names of all opened binary files */
static char* fbuff[_MAX_DBFIL]; /*!< In-core buffers for the binary files */
static size_t bptr[_MAX_DBFIL]; /*!< Indices to next available buffer location */
static size_t bsiz[_MAX_DBFIL]; /*!< The size of each in-core buffer */
static FT_FILE fdb[_MAX_DBFIL]; /*!< File pointers for each file */
static short doSwap[_MAX_DBFIL]; /*!< Byte-swapping flags for each file */
static const float Version = 1.0f; /*!< Internal file version tag */


SUBROUTINE(openbinarydb,OPENBINARYDB) (const char* filnam,
#ifdef _NCHAR_AFTER_CHARARG
                                       const int nchar,
#endif
                                       f90_int* filtyp, f90_int* ifile,
#ifdef _NCHAR_AFTER_CHARARG
                                       f90_int* status
#else
                                       f90_int* status, const int nchar
#endif
){
  int i, n;
  static char firstCall = 1;
  if (firstCall)
    for (i = 0; i < _MAX_DBFIL; i++)
    {
      fname[i] = fbuff[i] = NULL;
      bptr[i] = bsiz[i]  = 0;
      fdb[i] = 0;
      doSwap[0] = 0;
    }

  firstCall = 0;
  for (i = 0; i < _MAX_DBFIL; i++)
    if (!fdb[i] && (*filtyp > 0 || i > 0))
    {
      *ifile  = i;
      *status = 0;
      if (*filtyp > 0)
      {
        /* Get rid of trailing blanks and add a trailing 0-character instead */
        for (n = nchar; n > 0 && filnam[n-1] == ' '; n--);
        fname[i] = malloc(n+1);
        strncpy(fname[i],filnam,n);
        fname[i][n] = '\0';
      }
      switch (*filtyp)
      {
        case 0 : /* Temporary file */
#ifdef FT_USE_LOWLEVEL_IO
          fname[i] = tempnam("C:","fedem");
          fdb[i] = FT_open(fname[i],FT_tmp);
#else
          fdb[i] = tmpfile();
#endif
          break;
        case 1 : fdb[i] = FT_open(fname[i],FT_rb); break;
        case 2 : fdb[i] = FT_open(fname[i],FT_wb); break;
        case 3 : fdb[i] = FT_open(fname[i],FT_ab); break;
        case 4 : fdb[i] = FT_open(fname[i],FT_wbp); break;
        case 5 : fdb[i] = FT_open(fname[i],FT_rbp); break;
       default : *status = -abs(*filtyp);
      }
#ifdef FT_USE_LOWLEVEL_IO
      if (fdb[i] > 0)
        doSwap[i] = 0;
#else
      if (fdb[i])
        doSwap[i] = *filtyp > 0 ? 0 : -1;
#endif
      else if (*status < 0)
        fprintf(stderr,"openBinaryDB: Invalid file type %d\n",*filtyp);
      else
      {
        *status = -1;
        fprintf(stderr,"openBinaryDB: Error opening file of type %d\n",*filtyp);
        perror(fname[i] ? fname[i] : "temporary file");
      }
      if (*status < 0)
      {
        if (fname[i]) free(fname[i]);
        fname[i] = NULL;
        fdb[i] = 0;
      }
      return;
    }

  *status = -2;
  fprintf(stderr,"openBinaryDB: Too many binary files %d\n",i);
}


SUBROUTINE(setbufdb,SETBUFDB) (const f90_int* ifile, const f90_int* nbytes,
                               f90_int* status)
{
  *status = -1;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"setBufDB: Invalid file number %d\n",*ifile);
  else if (fbuff[*ifile])
    fprintf(stderr,"setBufDB: Binary file %d already has a buffer\n",*ifile);
  else if (*nbytes < 8)
    *status = 0; /* No buffer for this file */
  else if (!(fbuff[*ifile] = malloc(*nbytes)))
    fprintf(stderr,"setBufDB: Error allocating buffer for file %d\n",*ifile);
  else
  {
    bptr[*ifile] = 0;
    bsiz[*ifile] = *nbytes;
    *status = 0;
  }
}


/*!
  \brief Flushes the in-core buffer of the specified file to disk.
  \param[in] file Index of the file to flush
  \param[in] doFlush If true, perform a physical flush of the file
*/
static int flushBinaryDB (int file, int doFlush)
{
  size_t nWrote;
  if (fbuff[file] && bptr[file])
  {
    nWrote = FT_write(fbuff[file],1,bptr[file],fdb[file]);
    if (nWrote < bptr[file])
    {
      if (fname[file])
        fprintf(stderr,"flushBinaryDB: Error writing to file %s\n",fname[file]);
      else
        fprintf(stderr,"flushBinaryDB: Error writing to temporary file\n");
      fprintf(stderr,"               Wrote %llu of %llu bytes\n",
              (unsigned long long)nWrote,(unsigned long long)bptr[file]);
      return -2;
    }
    bptr[file] = 0;
  }
  else
    nWrote = 0;

  if (doFlush)
    if (FT_flush(fdb[file]))
    {
      if (fname[file])
        fprintf(stderr,"flushBinaryDB: Error flushing file %s\n",fname[file]);
      else
        fprintf(stderr,"flushBinaryDB: Error flushing temporary file\n");
      return -3;
    }

  return (nWrote > INT_MAX ? 1 : (int)nWrote);
}


SUBROUTINE(flushbinarydb,FLUSHBINARYDB) (const f90_int* ifile, f90_int* status)
{
  *status = -1;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"flushBinaryDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"flushBinaryDB: Binary file %d is not open\n",*ifile);
  else
    *status = flushBinaryDB (*ifile,1);
}


SUBROUTINE(closebinary_db,CLOSEBINARY_DB) (const f90_int* ifile,
                                           const f90_int* forceDelete,
                                           f90_int* status)
{
  *status = -1;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"closeBinaryDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"closeBinaryDB: Binary file %d is not open\n",*ifile);
  else if (flushBinaryDB(*ifile,0) >= 0 && FT_close(fdb[*ifile]))
    fprintf(stderr,"closeBinaryDB: Error closing binary file %d\n",*ifile);
  else
  {
    /* Delete a temporary file, and optionally also a regular file.
       If the deletetion of a regular file fails, this may be because
       it is simultaneously opened by another process (fedem_main).
       In that case, truncate the file instead to signal that it is void. */
    if (doSwap[*ifile] < 0 || *forceDelete)
      if (fname[*ifile] && unlink(fname[*ifile]) == -1 && *forceDelete)
        FT_close(FT_open(fname[*ifile],FT_wbp));

    if (fname[*ifile]) free(fname[*ifile]);
    if (fbuff[*ifile]) free(fbuff[*ifile]);
    fname[*ifile] = NULL;
    fbuff[*ifile] = NULL;
    bptr[*ifile] = 0;
    bsiz[*ifile] = 0;
    fdb[*ifile] = 0;
    *status = 0;
  }
}


SUBROUTINE(deletedb,DELETEDB) (const char* filnam, const int nchar)
{
  int n;
  char* fileName;
  for (n = nchar; n > 0 && filnam[n-1] == ' '; n--);
  fileName = malloc(n+1);
  strncpy(fileName,filnam,n);
  fileName[n] = '\0';
  unlink(fileName);
  free(fileName);
}


/*!
  \brief Closes all opened binary files.
  \details Used by signal handlers only.
*/
void closeAllBinaryDB ()
{
  int i;
  for (i = 0; i < _MAX_DBFIL; i++)
    if (fdb[i])
    {
      if (flushBinaryDB(i,0) >= 0 && FT_close(fdb[i]))
        fprintf(stderr,"closeAllBinaryDB: Error closing binary file %d\n",i);
      else
      {
        if (doSwap[i] < 0 && fname[i])
        {
          fprintf(stdout,"     Deleting temporary file %s\n",fname[i]);
          unlink(fname[i]);
        }
        if (fname[i]) free(fname[i]);
        if (fbuff[i]) free(fbuff[i]);
        fname[i] = NULL;
        fbuff[i] = NULL;
        bptr[i] = 0;
        bsiz[i] = 0;
        fdb[i] = 0;
      }
    }
}


/*!
  \brief Internal generic function to actually write binary data to file.

  \param[in] ifile File number in range [0,_MAX_DBFIL-1]
  \param[in] p     Pointer to the data to write
  \param[in] nSize Size of each data item
  \param[in] nData Number of data items to write
  \return 0 : Nothing is done (\a nData is zero)
  \return 1 : Number of bytes written &gt; \a INT_MAX
  \return &gt; 1 : Number of bytes written (&le; \a INT_MAX )
  \return -1 : Illegal file number or file not opened
  \return -2 : Error during write
*/
static int writeBinaryDB (int ifile, const void* p, size_t nSize, size_t nData)
{
  if (ifile < 0 || ifile >= _MAX_DBFIL)
    fprintf(stderr,"writeBinaryDB: Invalid file number %d\n",ifile);
  else if (!fdb[ifile])
    fprintf(stderr,"writeBinaryDB: Binary file %d is not open\n",ifile);
  else if (nData >= 1)
  {
    size_t nWrote;
    if (fbuff[ifile])
    {
      nWrote = nSize*nData;
      if (bptr[ifile] + nWrote > bsiz[ifile])
        if (flushBinaryDB(ifile,0) < 0) /* The file buffer is full */
          return -2;

      if (bptr[ifile] + nWrote <= bsiz[ifile])
      {
        /* Copy data to the file buffer */
        memcpy(fbuff[ifile]+bptr[ifile],p,nWrote);
        bptr[ifile] += nWrote;
        return (nWrote > INT_MAX ? 1 : (int)nWrote);
      }
    }

    /* Write data directly to file without buffering */
    nWrote = FT_write(p,nSize,nData,fdb[ifile]);
    if (nWrote < nData)
    {
      fprintf(stderr,"writeBinaryDB: Error writing to binary file %d\n",ifile);
      fprintf(stderr,"               Wrote %llu of %llu words\n",
              (unsigned long long)nWrote,(unsigned long long)nData);
      return -2;
    }

    nWrote *= nSize;
    return (nWrote > INT_MAX ? 1 : (int)nWrote);
  }
  else
    return 0;

  return -1;
}


/*!
  \brief Internal function to write double precision data to file.

  \param[in] ifile File number in range [0,_MAX_DBFIL-1]
  \param[in] data Pointer to the data to write
  \param[in] ndat Number of doubles to write

  \details This function is a wrapper for writeBinaryDB() used for
  large double arrays. The data is written in 2GB chunks on Windows
  to avoid run-time failure.
*/
static int writeDoubleDB (int ifile, const double* data, size_t ndat)
{
#if defined(win64)
  /* Note: There is a bug in the Microsoft C runtime library that results in
     fwrite hanging if the size of the data segment to write is 4GB or larger.
     See https://connect.microsoft.com/VisualStudio/feedback/details/755018/fwrite-hangs-with-large-size-count
     Thefore, we here write the data in 2GB chuncks to be on the safe side. */
  const size_t Giga = 2147483648; /* 2^31 = 2G */
#elif defined(aix64)
  /* Workaround for IBM bug */
  const size_t Giga = 1073741824; /* 2^30 = 1G */
#else
  const size_t Giga = 0;
#endif
  if (Giga > 0 && ndat > Giga/sizeof(double))
  {
    /* Write the data in GB segments */
    long long int nWrote, nTotal;
    size_t i, bufSiz;
    bufSiz = Giga/sizeof(double);
    nTotal = nWrote = 0;
    for (i = 0; i < ndat && nWrote >= 0; i += bufSiz)
    {
      if (i+bufSiz > ndat) bufSiz = ndat - i;
      nWrote = writeBinaryDB (ifile,(void*)(data+i),sizeof(double),bufSiz);
      if (nWrote < 0 || i == 0)
        nTotal = nWrote;
      else if (nTotal > 1)
        nTotal += nWrote;
    }
    return nTotal > INT_MAX ? 1 : (int)nTotal;
  }
  else
    return writeBinaryDB (ifile,(void*)data,sizeof(double),ndat);
}


SUBROUTINE(writeintdb,WRITEINTDB) (const f90_int* ifile, const f90_int* data,
                                   const f90_int* ndat, f90_int* status)
{
  *status = writeBinaryDB (*ifile,(void*)data,sizeof(f90_int),(size_t)*ndat);
}


SUBROUTINE(writefloatdb,WRITEFLOATDB) (const f90_int* ifile, const float* data,
                                       const f90_int* ndat, f90_int* status)
{
  *status = writeBinaryDB (*ifile,(void*)data,sizeof(float),(size_t)*ndat);
}


SUBROUTINE(writedoubled4,WRITEDOUBLED4) (const f90_int* ifile,
                                         const double* data,
                                         const f90_int* ndat, f90_int* status)
{
  *status = writeDoubleDB (*ifile,data,(size_t)*ndat);
}


SUBROUTINE(writedoubled8,WRITEDOUBLED8) (const f90_int* ifile,
                                         const double* data,
                                         const f90_int8* ndat, f90_int* status)
{
  *status = writeDoubleDB (*ifile,data,(size_t)*ndat);
}


SUBROUTINE(writechardb,WRITECHARDB) (const f90_int* ifile, char* data,
#ifdef _NCHAR_AFTER_CHARARG
                                     const int nchar, f90_int* status
#else
                                     f90_int* status, const int nchar
#endif
){
  if (data[nchar-1] == 13) data[nchar-1] = '\n';
  *status = writeBinaryDB (*ifile,(void*)data,sizeof(char),(size_t)nchar);
}


/*!
  \brief Performs byte swapping of an array in case of endian discrepancies.
*/
static int swapBytes (char* p, size_t m, size_t n)
{
  char* q;
  char* qend = p+m*n;

  switch (m)
  {
    case 2:
    {
      register char q0;
      for (q = p; q < qend; q += 2)
      {
        q0   = q[0];
        q[0] = q[1];
        q[1] = q0;
      }
      break;
    }
    case 4:
    {
      register char q0, q1;
      for (q = p; q < qend; q += 4)
      {
        q0   = q[0];
        q1   = q[1];
        q[0] = q[3];
        q[1] = q[2];
        q[2] = q1;
        q[3] = q0;
      }
      break;
    }
    case 8:
    {
      register char q0, q1, q2, q3;
      for (q = p; q < qend; q += 8)
      {
        q0   = q[0];
        q1   = q[1];
        q2   = q[2];
        q3   = q[3];
        q[0] = q[7];
        q[1] = q[6];
        q[2] = q[5];
        q[3] = q[4];
        q[4] = q3;
        q[5] = q2;
        q[6] = q1;
        q[7] = q0;
      }
      break;
    }
    default:
      fprintf(stderr,"swapBytes: Illegal word length %llu\n",
              (unsigned long long)m);
      return 0;
  }
  return (int)n;
}


/*!
  \brief Internal generic function to actually read binary data from file.

  \param[in] ifile File number in range &lt;-_MAX_DBFIL,_MAX_DBFIL&gt;. \n
                   If &lt; 0, rewind the file before reading it.
  \param[out] p    Pointer to memory segment in which to store the data read
  \param[in] nSize Size of each data item
  \param[in] nData Number of data items to read
  \param[in] silence If &gt; 0, suppress error message on read failure
  \return 0 : Nothing is done (\a nData is zero)
  \return &gt; 0 : Number of bytes read
  \return -1 : Illegal file number or file not opened
  \return -2 : Error during read
  \return -3 : Invalid record length
  \return -99 : End-of-file detected (no data read)
*/
static int readBinaryDB (int ifile, void* p, size_t nSize,
                         size_t nData, char silence)
{
  if (!(ifile < _MAX_DBFIL && -ifile < _MAX_DBFIL))
    fprintf(stderr,"readBinaryDB: Invalid file number %d\n",ifile);
  else if (!fdb[ifile >= 0 ? ifile : -ifile])
    fprintf(stderr,"readBinaryDB: Binary file %d is not open\n",ifile);
  else if (nData >= 1)
  {
    size_t nRead;
    if (ifile < 0) FT_seek(fdb[-ifile],(FT_int)0,SEEK_SET);
    nRead = FT_read(p,nSize,nData,fdb[ifile >= 0 ? ifile : -ifile]);
    if (nRead < nData)
    {
      const char* filnam = fname[ifile >= 0 ? ifile : -ifile];
      if (silence)
        return nRead == 0 ? -99 : -2;
      else if (filnam)
        fprintf(stderr,"readBinaryDB: Error reading from file %s\n",filnam);
      else
        fprintf(stderr,"readBinaryDB: Error reading from file %d\n",ifile);
      fprintf(stderr,"              Read %llu of %llu words\n",
              (unsigned long long)nRead,(unsigned long long)nData);
      return -2;
    }
    else if (doSwap[ifile >= 0 ? ifile : -ifile] > 0 && nSize > 1)
      if (!swapBytes((char*)p,nSize,nData))
        return -3;

    nRead *= nSize;
    return (nRead > INT_MAX ? 1 : (int)nRead);
  }
  else
    return 0;

  return -1;
}


SUBROUTINE(readintdb,READINTDB) (const f90_int* ifile, f90_int* data,
                                 const f90_int* ndat, f90_int* status)
{
  if (*ndat == -1) /* checking for more data, no message on read failure */
    *status = readBinaryDB (*ifile,(void*)data,sizeof(f90_int),1,1);
  else
    *status = readBinaryDB (*ifile,(void*)data,sizeof(f90_int),(size_t)*ndat,0);
}


SUBROUTINE(readfloatdb,READFLOATDB) (const f90_int* ifile, float* data,
                                     const f90_int* ndat, f90_int* status)
{
  *status = readBinaryDB (*ifile,(void*)data,sizeof(float),(size_t)*ndat,0);
}


SUBROUTINE(readdoubled4,READDOUBLED4) (const f90_int* ifile, double* data,
                                       const f90_int* ndat, f90_int* status)
{
  *status = readBinaryDB (*ifile,(void*)data,sizeof(double),(size_t)*ndat,0);
}

SUBROUTINE(readdoubled8,READDOUBLED8) (const f90_int* ifile, double* data,
                                       const f90_int8* ndat, f90_int* status)
{
  *status = readBinaryDB (*ifile,(void*)data,sizeof(double),(size_t)*ndat,0);
}


SUBROUTINE(readchardb,READCHARDB) (const f90_int* ifile, char* data,
#ifdef _NCHAR_AFTER_CHARARG
                                   const int nchar, const f90_int* ndat,
                                   f90_int* status
#else
                                   const f90_int* ndat, f90_int* status,
                                   const int nchar
#endif
){
  size_t nbytes = (size_t)(*ndat < nchar ? *ndat : nchar);
  *status = readBinaryDB (*ifile,(void*)data,sizeof(char),nbytes,0);
  if (*ndat < nchar && *status >= 0) memset(data+(*ndat),' ',nchar-(*ndat));
}


SUBROUTINE(writetagdb,WRITETAGDB) (const f90_int* ifile, const char* tag,
#ifdef _NCHAR_AFTER_CHARARG
                                   const int nchar, const f90_int* cs,
                                   f90_int* status
#else
                                   const f90_int* cs, f90_int* status,
                                   const int nchar
#endif
){
  FT_int start;
  char buf[16];
  size_t ncharB;
  sprintf(buf,";%.1f;\n",Version);
  ncharB = strlen(buf);

  *status = -99;
  start = *ifile < 0 || *ifile >= _MAX_DBFIL ? (FT_int)0 : FT_tell(fdb[*ifile]);
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"writeTagDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"writeTagDB: Binary file %d is not open\n",*ifile);
  else if (FFa_writeTag(fdb[*ifile],tag,nchar,(unsigned int)*cs) < 0)
    fprintf(stderr,"writeTagDB: Error writing file tag to %s\n",fname[*ifile]);
  else if (FT_write(buf,1,ncharB,fdb[*ifile]) < ncharB)
    fprintf(stderr,"writeTagDB: Error writing version # to %s\n",fname[*ifile]);
  else
    *status = (f90_int)(FT_tell(fdb[*ifile]) - start);
}


SUBROUTINE(readtagdb,READTAGDB) (const f90_int* ifile, char* tag,
#ifdef _NCHAR_AFTER_CHARARG
                                 const int nchar, f90_int* cs, f90_int* status
#else
                                 f90_int* cs, f90_int* status, const int nchar
#endif
){
  FT_int start;
  size_t len_Tag;
  int nread, nval;
  float currVer, iptr;
  char buf[16];

  start = 0;
  *status = -99;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"readTagDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"readTagDB: Binary file %d is not open\n",*ifile);
  else if (doSwap[*ifile] < 0)
    fprintf(stderr,"readTagDB: File %s is a temporary file\n",fname[*ifile]);
  else
  {
    start = FT_tell(fdb[*ifile]);
    *status = FFa_readTag(fdb[*ifile],tag,nchar,(unsigned int*)cs);
  }

  if (*status == 0)
  {
    fprintf(stderr,"readTagDB: File %s is not a binary file\n",fname[*ifile]);
    *status = -9;
    return;
  }
  else if (*status < 0)
    return;

  /* Check for byte swapping */
  doSwap[*ifile] = (*status != FFa_endian());
  if (doSwap[*ifile]) swapBytes((char*)cs,sizeof(unsigned int),1);

  /* Fill tag variable with trailing blanks if necessary */
  len_Tag = strlen(tag);
  if ((int)len_Tag < nchar) memset(tag+len_Tag,' ',nchar-len_Tag);

  /* Check file version number */
  *status = -98;
  if (!FT_gets(buf,16,fdb[*ifile]))
    fprintf(stderr,"readTagDB: Error reading file version number\n");
  else if ((nread = sscanf(buf,";%f;%d;",&currVer,&nval)) < 1)
    fprintf(stderr,"readTagDB: Erroneous file version field: %s",buf);
  else if (currVer != Version && nread < 2)
    fprintf(stderr,"readTagDB: Wrong file version: %f\n",currVer);
  else if (nread == 2 && currVer > 2.0f && nval > 0)
    /* Return the value stored after the version tag + the minor version */
    *status = 10*nval + (int)(10.0f*modff(currVer,&iptr) + 0.5f);
  else
    *status = (f90_int)(FT_tell(fdb[*ifile]) - start);
}


SUBROUTINE(setpositiondb,SETPOSITIONDB) (const f90_int* ifile,
                                         const f90_int8* pos, f90_int* status)
{
  *status = -1;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"setPositionDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"setPositionDB: Binary file %d is not open\n",*ifile);
  else if (FT_seek(fdb[*ifile],(FT_int)*pos,SEEK_SET) < 0)
    fprintf(stderr,"setPositionDB: Improper seek on file %s\n",fname[*ifile]);
  else
    *status = 0;
}


SUBROUTINE(getpositiondb,GETPOSITIONDB) (const f90_int* ifile, f90_int8* pos)
{
  FT_int newPos;
  newPos = (FT_int)0;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"getPositionDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"getPositionDB: Binary file %d is not open\n",*ifile);
  else
    newPos = FT_tell(fdb[*ifile]);
#if defined(FT_USE_LOWLEVEL_IO) && defined(FT_HAS_INT4_ONLY)
  if (newPos > INT_MAX)
  {
    fprintf(stderr,"getPositionDB: Binary file too large %s\n",fname[*ifile]);
    fprintf(stderr,"               Cannot do random access over 2GB files.\n");
    *pos = 0;
    return;
  }
#endif
  *pos = newPos;
}


SUBROUTINE(putchardb,PUTCHARDB) (const f90_int* ifile, const char* tag,
#ifdef _NCHAR_AFTER_CHARARG
                                 const int ncharT, const char* data,
                                 const int ncharD, f90_int* status
#else
                                 const char* data, f90_int* status,
                                 const int ncharT, const int ncharD
#endif
){
  char* buf;
  int i;

  *status = -1;
  if (*ifile < 0 || *ifile >= _MAX_DBFIL)
    fprintf(stderr,"putCharDB: Invalid file number %d\n",*ifile);
  else if (!fdb[*ifile])
    fprintf(stderr,"putCharDB: Binary file %d is not open\n",*ifile);
  else
    *status = 0;
  if (*status) return;

  FT_seek(fdb[*ifile],(FT_int)0,SEEK_SET);

  /* Fill buf with the first ncharT bytes from the file */
  buf = malloc(ncharT+1);
  memset(buf,0,ncharT+1);
  for (i = 0; i < ncharT && !FT_eof(fdb[*ifile]); i++)
    buf[i] = FT_getc(fdb[*ifile]);

  /* Continue reading one-by-one byte from the file adding it to the end of buf
     while shifting the ncharT-1 earlier bytes one step forward */
  while (!FT_eof(fdb[*ifile]))
  {
    if (!strncmp(buf,tag,ncharT))
    {
      /* Yes, current contents of buf matches our tag.
         Now overwrite the data following immediately after the tag */
      FT_seek(fdb[*ifile],(FT_int)0,SEEK_CUR);
      *status = (f90_int)FT_write(data,1,ncharD,fdb[*ifile]) - (f90_int)ncharD;
      if (*status < 0)
        fprintf(stderr,"putCharDB: Error writing to file %s\n",fname[*ifile]);
      free(buf);
      return;
    }

    /* Shift buf one byte forward and read the next byte from file */
    for (i = 1; i < ncharT; i++) buf[i-1] = buf[i];
    buf[ncharT-1] = FT_getc(fdb[*ifile]);
  }

  strncpy(buf,tag,ncharT);
  fprintf(stderr,"putCharDB: End-of-file during read of %s\n",fname[*ifile]);
  fprintf(stderr,"           Could not find tag %s\n",buf);
  free(buf);
  *status = 1;
}


SUBROUTINE(copybinarydb,COPYBINARYDB) (const f90_int* tofile,
                                       const f90_int* fromfile,
                                       const f90_int* nBytes,
                                       f90_int* status)
{
  void* buffer;

  *status = 0;
  if (*nBytes < 1) return;

  buffer = malloc(*nBytes);
  *status = readBinaryDB (-(*fromfile),buffer,1,(size_t)*nBytes,0);
  if (*status > 0)
    *status = writeBinaryDB (*tofile,buffer,1,(size_t)*status);

  free(buffer);
}
