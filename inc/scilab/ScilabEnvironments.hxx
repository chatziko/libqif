/*
 * Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
 * Copyright (C) 2012 - Scilab Enterprises - Calixte DENIZET
 *
 * This file must be used under the terms of the CeCILL.
 * This source file is licensed as described in the file COPYING, which
 * you should have received as part of this distribution.  The terms
 * are also available at
 * http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
 *
 */

#ifndef __SCILABENVIRONMENTS_HXX__
#define __SCILABENVIRONMENTS_HXX__

#include <vector>

#include "ScilabAbstractEnvironmentException.hxx"
#include "ScilabAbstractEnvironment.hxx"
#include "dynlib_external_objects_scilab.h"

extern "C" {
#include "localization.h"
}

namespace org_modules_external_objects
{

class EXTERNAL_OBJECTS_SCILAB_IMPEXP ScilabEnvironments
{
    static std::vector<ScilabAbstractEnvironment*> environments;

public:

    static int registerScilabEnvironment(ScilabAbstractEnvironment * env);

    static void unregisterScilabEnvironment(const int id);

    static ScilabAbstractEnvironment & getEnvironment(const int id);
};

}

#endif // __SCILABENVIRONMENTS_HXX__
