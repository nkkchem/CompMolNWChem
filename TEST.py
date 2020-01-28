import errno
import json
import multiprocessing
import os
import re
import subprocess
import sys
import time
import traceback
import uuid
import zipfile

os.system('ls')
result_directory = './data'

output_files = list()
result_file = os.path.join(result_directory, 'TopHat2_result.zip')

print(result_file)

result_dirs = os.listdir(result_directory)

print(result_dirs)

tophat2_result_dir_name = filter(re.compile('tophat2_result_*').match, result_dirs)[0]
tophat2_result_dir = os.path.join(result_directory, tophat2_result_dir_name)


