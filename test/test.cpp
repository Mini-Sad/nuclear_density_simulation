/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/TestMain.h>
#include <cxxtest/ErrorPrinter.h>

int main( int argc, char *argv[] ) {
 int status;
    CxxTest::ErrorPrinter tmp;
    CxxTest::RealWorldDescription::_worldName = "cxxtest";
    status = CxxTest::Main< CxxTest::ErrorPrinter >( tmp, argc, argv );
    return status;
}
bool suite_TestPoly_init = false;
#include "test_poly.h"

static TestPoly suite_TestPoly;

static CxxTest::List Tests_TestPoly = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestPoly( "test_poly.h", 4, "TestPoly", suite_TestPoly, Tests_TestPoly );

static class TestDescription_suite_TestPoly_testPoly : public CxxTest::RealTestDescription {
public:
 TestDescription_suite_TestPoly_testPoly() : CxxTest::RealTestDescription( Tests_TestPoly, suiteDescription_TestPoly, 8, "testPoly" ) {}
 void runTest() { suite_TestPoly.testPoly(); }
} testDescription_suite_TestPoly_testPoly;

#include <cxxtest/Root.cpp>
const char* CxxTest::RealWorldDescription::_worldName = "cxxtest";
