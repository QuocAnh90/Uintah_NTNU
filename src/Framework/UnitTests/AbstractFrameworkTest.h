/*
  For more information, please see: http://software.sci.utah.edu

  The MIT License

  Copyright (c) 2006 Scientific Computing and Imaging Institute,
  University of Utah.

  License for the specific language governing rights and limitations under
  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included
  in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.
*/

#ifndef Framework_UnitTests_AbstractFrameworkTest_h
#define Framework_UnitTests_AbstractFrameworkTest_h

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Framework/SCIRunFramework.h>
#include <Core/CCA/spec/cca_sidl.h>

class AbstractFrameworkTest : public CppUnit::TestFixture {
public:
  AbstractFrameworkTest();
  virtual ~AbstractFrameworkTest();

  // Set up context before running a test.
  virtual void setUp();

  // Clean up after the test run.
  virtual void tearDown();

  void testInstantiate();
  void testGetServices();

private:
  CPPUNIT_TEST_SUITE( AbstractFrameworkTest );
  CPPUNIT_TEST( testInstantiate );
  CPPUNIT_TEST( testGetServices );
  CPPUNIT_TEST_SUITE_END();

//   sci::cca::AbstractFramework::pointer abstractFramework;
//   SCIRun::SCIRunFramework *fwk;
};

#endif
