/*
 * Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
 * Copyright (C) INRIA
 *
 * This file must be used under the terms of the CeCILL.
 * This source file is licensed as described in the file COPYING, which
 * you should have received as part of this distribution.  The terms
 * are also available at
 * http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
 *
 */
#ifndef SCI_VERSION_H
#define SCI_VERSION_H

#define SCI_VERSION_MAJOR 5
#define SCI_VERSION_MINOR 4
#define SCI_VERSION_MAINTENANCE 1
#define SCI_VERSION_STRING "scilab-5.4.1"
/* SCI_VERSION_REVISION --> hash key commit */
#define SCI_VERSION_SVN_URL "svn://svn.scilab.org/scilab/branches/5.4/"/
#define SCI_VERSION_GIT_URL "git://git.scilab.org/scilab/"
#define SCI_VERSION_REVISION origin/5.4
#define SCI_VERSION_TIMESTAMP 1364571296

void disp_scilab_version(void);

/* for compatibility */
/* Deprecated */
#define SCI_VERSION SCI_VERSION_STRING
#define DEFAULT_SCI_VERSION_MESSAGE "scilab-5.4.1 (INRIA,ENPC)"


#endif
/*--------------------------------------------------------------------------*/

