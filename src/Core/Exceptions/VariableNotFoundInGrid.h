/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */


/*
 *  VariableNotFoundInGrid.h: 
 *
 *  Written by:
 *   James Bigler
 *   Department of Computer Science
 *   University of Utah
 *   May 2001
 *
 */

#ifndef UINTAH_EXCEPTIONS_VARIABLENOTFOUNDINGRID_H
#define UINTAH_EXCEPTIONS_VARIABLENOTFOUNDINGRID_H

#include <Core/Exceptions/Exception.h>
#include <Core/Geometry/IntVector.h>
#include <string>

namespace Uintah {
  
  class VariableNotFoundInGrid : public Uintah::Exception {
  public:
    VariableNotFoundInGrid(const std::string& varname, long particleID,
                           int matlIndex, const std::string& extramsg, 
                           const char* file, int line);
    VariableNotFoundInGrid(const std::string& varname, IntVector loc,
                           int matlIndex, const std::string& extramsg,
                           const char* file, int line);
    VariableNotFoundInGrid(const std::string& varname,
                           const std::string& extramsg,
                           const char* file, int line);
    VariableNotFoundInGrid(const VariableNotFoundInGrid&);
    virtual ~VariableNotFoundInGrid();
    virtual const char* message() const;
    virtual const char* type() const;
  protected:
  private:
    std::string d_msg;
    VariableNotFoundInGrid& operator=(const VariableNotFoundInGrid&);
  };
  
} // End namespace Uintah

#endif


