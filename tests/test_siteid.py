import unittest
from test_tools import cd, run_integration_test

class TestSiteID( unittest.TestCase ):

    def test_siteid_runs_with_old_data( self ):
        data_path = 'siteid_tests/garnet_tet_and_oct'
        executable = 'siteid'
        run_integration_test( data_path, executable )
    
    def test_siteid_finds_spherical_site( seld ):
        data_path = 'siteid_tests/find_single_spherical_site'                
        executable = 'siteid'
        run_integration_test( data_path, executable )
        
if __name__ == '__main__':
    unittest.main()
