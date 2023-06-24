from getGC import *

def test_empty_string():
	s = ""
	assert getGC([s]) == (0.0,"",0)

def test_nonGC_string():
	s = "abtbabfabfabagbab"
	assert getGC([s]) == (0.0, "abtbabfabfabagbab", 0)