import unittest
from test_tools import cd, run_integration_test

class TestSiteID( unittest.TestCase ):

    def test_polyhedra_finds_tetrahedra( self ):
        data_path = 'find_polyhedra_tests/find_tetrahedra'
        executable = 'find_polyhedra'
        run_integration_test( data_path, executable )
    
if __name__ == '__main__':
    unittest.main()
