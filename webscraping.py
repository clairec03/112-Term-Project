from bs4 import BeautifulSoup
import requests
input = "mongoose"
url = f"https://www.google.com/search?q={input}"

page = requests.get(url).text
soup = BeautifulSoup(page, "html.parser")
imageClass = soup.find('img', attrs={"class": "rISBZc M4dUYb"})
print(imageClass)
try: text = imageClass.get_text()
except: 
    print("Image not found")
print(text)