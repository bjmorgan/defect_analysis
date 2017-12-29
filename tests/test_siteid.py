import unittest
from test_tools import cd, run_integration_test

class TestSiteID( unittest.TestCase ):

    def test_siteid_runs_with_old_data( self ):
        data_path = 'test_siteid'
        executable = 'siteid'
        run_integration_test( data_path, executable )
    
    def test_siteid_finds_spherical_site( seld ):
        data_path = 'test_sphereid'                
        executable = 'siteid'
        run_integration_test( data_path, executable )
        
if __name__ == '__main__':
    unittest.main()
