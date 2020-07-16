#Connect to RSS feeds to create and update component databases
import feedparser
import os
import shutil

databaseRssUrls = ["http://www.rc-international.info/rss/","http://www.rc-international.info/rss/arrivals/","https://hobbyking.com/en_us/rss/catalog/new/store_id/1/"]

for url in databaseRssUrls:
    rssParse = feedparser.parse(url)
    print("\n\n")
    print(rssParse["feed"]["title"])
    print("\n")

    for entry in rssParse["entries"]:
        
        entryTitle = entry["title"]
        print(entryTitle)
        entryTitle = entryTitle.lower()
        entryType = None
        
        if ("battery" in entryTitle and not "charger" in entryTitle) or "mah" in entryTitle:
            entryType = "Battery"
            print(entryType)
            continue
        
        if "esc" in entryTitle and not "fan" in entryTitle:
            entryType = "ESC"
            print(entryType)
            continue