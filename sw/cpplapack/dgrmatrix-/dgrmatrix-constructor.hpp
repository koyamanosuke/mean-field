//=============================================================================
/*! dgrmatrix constructor without arguments */
inline dgrmatrix::dgrmatrix()
{CPPL_VERBOSE_REPORT;
  m =0;
  n =0;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgrmatrix copy constructor */
inline dgrmatrix::dgrmatrix(const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  copy(mat);
}

//=============================================================================
/*! dgrmatrix constructor with filename */
inline dgrmatrix::dgrmatrix(const char* filename)
{CPPL_VERBOSE_REPORT;
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgrmatrix destructor */
inline dgrmatrix::~dgrmatrix()
{CPPL_VERBOSE_REPORT;
  a.clear();
  ia.clear();
  ja.clear();
}
