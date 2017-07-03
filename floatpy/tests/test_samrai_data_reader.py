import os
import unittest

from floatpy.readers import samrai_reader

class TestSamraiDataReader(unittest.TestCase):

    def setUp(self):
        self.directory_name = os.path.join(os.path.dirname(__file__), 'test_samrai_data_reader_data')
        data_reader = samrai_reader.SamraiDataReader(self.directory_name)

if __name__ == '__main__':
    unittest.main()
