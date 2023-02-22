
import traceback
import logging

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def print_test_summary(file_name, test_dict):

	print("\n%s TEST SUMMARY:"%file_name)
	
	for test in test_dict:
		if test_dict[test] is None:
			result = f"{bcolors.OKGREEN}PASSED{bcolors.ENDC}"
		else: 
			result = f"{bcolors.FAIL}FAILED{bcolors.ENDC}"
		print(test, ":", result)

		if test_dict[test] is not None:
			print(test_dict[test])

	print()

