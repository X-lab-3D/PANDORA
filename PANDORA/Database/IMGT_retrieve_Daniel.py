import mechanize

BASE_URL = 'https://www.imgt.org/3Dstructure-DB/'

def createNewBrowser():
	# Create browser
	browser = mechanize.Browser()
	# Makes site know I am not a robot ;)
	browser.set_handle_robots(False) 
	return browser
	
def retrieveHtml(browser):
	# Open site, select fitness and retrieve the html code
	browser.select_form(nr=0)
	browser.form.find_control(name="type-entry").value = ["PDB"]
	browser.form.find_control(name="ReceptorType").value = ["peptide/MH2"]
	browser.submit()
	html = browser.response().read().decode("utf-8")
	return html

if __name__ == "__main__":
	browser = createNewBrowser()
	browser.open(BASE_URL)
	# get the html
	rawhtml = retrieveHtml(browser)
	# save the html
	writer = open('html.html', 'w')
	writer.write(rawhtml)
	writer.close()
