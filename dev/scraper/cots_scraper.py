import requests
import json
from selenium import webdriver

if __name__=="__main__":

    # Hobby King
    #url = "https://hobbyking.com/en_us/batteries-chargers/batteries.html#q=&idx=hbk_live_magento_en_us_products&dFR[warehouses][0]=USA&dFR[warehouses][1]=Global&dFR[warehouses_stock_data][0]=USA|1&dFR[warehouses_stock_data][1]=USA|2&dFR[warehouses_stock_data][2]=USA|3&dFR[warehouses_stock_data][3]=Global|1&dFR[warehouses_stock_data][4]=Global|2&dFR[warehouses_stock_data][5]=Global|3&hFR[categories.level0][0]=Batteries %2F Chargers %2F%2F%2F Batteries&is_v=1"
    url = "https://hobbyking.com/en_us/batteries-chargers/batteries.html"

    # Initialize driver
    path = "/usr/bin/firefox"
    driver = webdriver.Firefox(executable_path=path)

    # Get webpage
    driver.get(url)