import tempfile
import os
import shutil
import subprocess
import numpy as np

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def run_integration_test( data_path, executable ):
    test_dir = tempfile.mkdtemp()
    src = './{}/inputs'.format( data_path )
    input_files = os.listdir( src )
    for f in input_files:
       full_file_name = os.path.join(src, f)
       if (os.path.isfile(full_file_name)):
           shutil.copy( full_file_name, test_dir )
    shutil.copy( '../bin/{}'.format( executable ), test_dir )
    with cd( test_dir ):
        p = subprocess.run( [ './{}'.format( executable ) ],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE )
    output_src = './{}/expected_outputs'.format( data_path )
    output_files = os.listdir( output_src )
    for of in output_files:
        expected_data = np.loadtxt( os.path.join( output_src, of ) )
        calculated_data = np.loadtxt( os.path.join( test_dir, of ) )
        np.testing.assert_array_equal( expected_data, calculated_data )
