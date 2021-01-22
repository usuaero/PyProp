# Define here the models for your scraped items
#
# See documentation in:
# https://docs.scrapy.org/en/latest/topics/items.html

import scrapy


class TowerHobbiesItem(scrapy.Item):
    # define the fields for your item here like:
    # name = scrapy.Field()
    Title = scrapy.Field()
    Manufacturer = scrapy.Field()
    Bat_Capacity = scrapy.Field()
    Bat_DischargeRate = scrapy.Field()
    Bat_Weight = scrapy.Field()
    Bat_Voltage = scrapy.Field()
    Bat_Config = scrapy.Field()
    Bat_Resistance = scrapy.Field()
    Bat_Chemistry = scrapy.Field()
    Mot_Kv = scrapy.Field()
    Mot_Weight = scrapy.Field()
    Mot_NoLoadCurrent = scrapy.Field()
    Mot_Resistance = scrapy.Field()
    Mot_MaxCurrent = scrapy.Field()
    Mot_GearRatio = scrapy.Field()
    ESC_MaxCurrent = scrapy.Field()
    ESC_Resistance = scrapy.Field()
    ESC_Weight = scrapy.Field()
    
    pass
