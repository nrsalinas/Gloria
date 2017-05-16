import unittest
import gloria_test

def my_module_suite():
	loader = unittest.TestLoader()
	suite = loader.loadTestsFromModule(gloria_test)
	return suite
