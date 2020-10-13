# imports the libraries "requests" used for webscraping and "html" from "lxml" to read html sources
from lxml import html
import requests

# webscraper function inputs a website and then pulls the information and prints it neatly
def webscraper(website):
    # defines variable to use a certain website url
    page = requests.get(website)
    tree = html.fromstring(page.content)

    # variable for defining the 'full XML' path for the desired info
    arrayTest = "/html/body/div[1]/div/div[3]/div[1]/form/div/div[3]/div[1]/div[2]/div/div/div/text()"

    # prints the title for the desired information e.g "SPECIFICATIONS" from the actual website
    titleStr = tree.xpath('/html/body/div[1]/div/div[3]/div[1]/form/div/div[3]/div[1]/div[2]/div/div/div/b[4]')[0].text
    print("\n", titleStr)

    # creates a variable for iterating through the information needed
    intList = 6

    # The info needed to pull the information from the website is a list element in the form of a string in the
    # variable "arrayTest" above. e.g. .../div/div/div/text()[0], .../div/div/div/text()[1], .../div/div/div/text()[2]
    # this for loop below increments the amount of lines I want to pull and print from the html website directories
    # provided from the variable arrayTest.
    for i in range(0, 4):
        i += intList

        #creates a string that looks like a list e.g. "[0]", "[2]"
        intToStr = "[" + str(i) + "]"

        # adds both the variable arrayTest to intToStr. e.g. .../div/text() + [0] = .../div/text()[0].
        # the format .../div/text()[0] is then readable by the .xpath() method in python
        motorSpecs = tree.xpath(arrayTest + intToStr)

        # prints the readable variable "motorSpecs" and formats it nicely using .strip() methods.
        print(str(motorSpecs).strip('[]').strip('\'\'').strip())


if __name__ == '__main__':
    # input website for pulling information needed through "webscraper()" method
    websiteHtml = 'https://vortexhobbies.com/traxxas-motor-titan-12t-550-turns-tra3785-p-607.html'
    webscraper(websiteHtml)


