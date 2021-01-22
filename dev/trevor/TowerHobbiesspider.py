# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 23:51:23 2020

@author: Trevor Coombs
"""

import scrapy
from bs4 import BeautifulSoup
import re
from TowerHobbies.items import TowerHobbiesItem


class TowerHobbiesSpider(scrapy.Spider):
    name = 'TowerHobbiesspider'
    allowed_domains = ['towerhobbies.com']
    start_urls = ['https://www.towerhobbies.com/']

    def __init__(self):
        self.declare_xpath()

        #All the XPaths the spider needs to know go here
    def declare_xpath(self):
        self.getAllCategoriesXpath = "//*[@id='sg-navbar-collapse']/div/div/nav/div[2]/ul/li[1]/a"
        self.getAllSubCategoriesXpath = "//*[@id='sg-navbar-collapse']/div/div/nav/div[2]/ul/li[1]/div/div/div/div[2]/ul/li/a"
        self.getAllSubCategoriesTwoXpath = "//*[@id='sg-navbar-collapse']/div/div/nav/div[2]/ul/li[1]/div/div/div/div[2]/ul/li/div/ul/li[4]/a"
        self.getAllSubCategoriesThreeXpath = "//*[@id='sg-navbar-collapse']/div/div/nav/div[2]/ul/li[1]/div/div/div/div[2]/ul/li/div/ul/li[4]/ul/li[1]/a"
        self.getAllItemsXpath = "//*[@id='product-search-results']/div/div[2]/div/div[2]/div[2]/div[1]/div/div[1]/div[2]/div[3]/a"
        self.TitleXpath  = "//*[@id='maincontent']/div[1]/div/div[1]/div/div[2]/div/div/h1/text()"
        self.ManufacturerXpath = ""
        self.BatteryXpath = ""
        self.MotorXpath = ""
        self.ESCXpath = ""
      #  self.SpecsXpath = ""

    def parse(self, response):
        for href in response.xpath(self.getAllCategoriesXpath):
            url = response.urljoin(href.extract())
            yield scrapy.Request(url=url,callback=self.parse_category)
 
    def parse_category(self,response):
        for href in response.xpath(self.getAllSubCategoriesXpath):
            url = response.urljoin(href.extract())
            yield scrapy.Request(url,callback=self.parse_subcategory)
    
    def parse_categoryTwo(self,response):
        for href in response.xpath(self.getAllSubCategoriesTwoXpath):
            url = response.urljoin(href.extract())
            yield scrapy.Request(url,callback=self.parse_subcategory)
            
    def parse_categoryThree(self,response):
        for href in response.xpath(self.getAllSubCategoriesTwoXpath):
            url = response.urljoin(href.extract())
            yield scrapy.Request(url,callback=self.parse_subcategory)

    def parse_subcategory(self,response):
        for href in response.xpath(self.getAllItemsXpath):
            url = response.urljoin(href.extract())
            yield scrapy.Request(url,callback=self.parse_main_item)

    
    def parse_main_item(self,response):
        item = TowerHobbiesItem()
 
        Title = response.xpath(self.TitleXpath).extract()
        Title = self.cleanText(self.parseText(self.listToStr(Title)))
 
        Manufacturer = response.xpath(self.ManufacturerXpath).extract()
        Manufacturer = self.cleanText(self.parseText(self.listToStr(Manufacturer)))
 
        Battery = response.xpath(self.BatteryXpath).extract()
        Battery = self.cleanText(self.parseText(self.listToStr(Battery)))
 
        Motor = response.xpath(self.MotorXpath).extract()
        Motor = self.cleanText(self.parseText(self.listToStr(Motor)))
 
        ESC = response.xpath(self.ESCXpath).extract()
        ESC = self.cleanText(self.parseText(self.listToStr(ESC)))

        #Specs = response.xpath(self.SpecsXpath).extract()
        #Specs = self.cleanText(self.parseText(Specs))

        #Put each element into its item attribute.
        item['Title']          = Title
        item['Manufacturer']   = Manufacturer
        item['Battery']         = Battery
        item['Motor']           = Motor
        item['ESC']             = ESC
        #item['Specs']           = Specs
        return item
 
    #Methods to clean and format text to make it easier to work with later
    def listToStr(self,MyList):
        dumm = ""
        MyList = [i.encode('utf-8') for i in MyList]
        for i in MyList:dumm = "{0}{1}".format(dumm,i)
        return dumm
 
    def parseText(self, str):
        soup = BeautifulSoup(str, 'html.parser')
        return re.sub(" +|\n|\r|\t|\0|\x0b|\xa0",' ',soup.get_text()).strip()
 
    def cleanText(self,text):
        soup = BeautifulSoup(text,'html.parser')
        text = soup.get_text()
        text = re.sub("( +|\n|\r|\t|\0|\x0b|\xa0|\xbb|\xab)+",' ',text).strip()
        return text