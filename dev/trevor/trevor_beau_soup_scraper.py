from urllib.request import urlopen as uReq
from bs4 import BeautifulSoup as soup

tower_hobbies_batteries = 'https://www.towerhobbies.com/search?q=Battery&search-button=&lang=default'

#opening up connection, grabbing the page
uClient = uReq(tower_hobbies_batteries)
#offloads content into a variable
batteries_html = uClient.read()
#closes the client
uClient.close()



#Does html parsing
page_soup = soup(batteries_html, "html.parser")

#Grabs each product
tiles = page_soup.findAll("div",{"class":"product-tile"})

for tile in tiles:
    print(tile)

# There are a lot of links you have to go through to get to the actual product pages.