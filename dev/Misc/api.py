#Connect to APIs to create and update component databases
import os
import shutil
import requests

response = requests.get("http://www.hobbyking.com/hobbyking_api")
print(response.status_code)