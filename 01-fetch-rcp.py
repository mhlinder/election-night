
import pdb
import re

import bs4
import pandas as pd
import requests

url_host = 'https://www.realclearpolitics.com'
url_base = '{}/epolls/latest_polls/{{}}'.format(url_host)

outdir = 'data'

def get_soup(url):
    r = requests.get(url)
    soup = bs4.BeautifulSoup(r.text)
    return soup

def get_races(race_types = ['senate', 'house', 'governor']):
    races = {}
    for race in race_types:
        soup = get_soup(url_base.format(race))
        links = []
        for td in soup.find_all('td', attrs={'class': 'lp-race'}):
            a = td.find('a')
            link = '{}{}'.format(url_host, a['href'])
            links.append(link)
        links = set(links)
        races[race] = links
    return races

def get_race_data(url, start_year = 2018):
    soup = get_soup(url)
    tables = soup.find_all('table', attrs={'class': ['data', 'large']})

    data = []
    for table in tables:
        data.append(table.find_all('tr'))

    rows = []
    for table in data:
        last_month = None
        header = []
        k = 0
        for row in table:
            if k == 0:
                year = start_year
            else:
                last_month = this_month
                
            if 'class' in row.attrs.keys() and 'header' in row['class']:
                for th in row.find_all('th'):
                    header.append(th.text)
                continue
            
            newdata = []
            
            tds = row.find_all('td')
            if tds[0].text =='RCP Average':
                continue

            newdata.append(tds[0].find(attrs={'class': 'normal_pollster_name'}).text)
            for td in tds[1:]:
                newdata.append(td.text)
            newrow = dict(zip(header, newdata))

            this_month = int(newrow['Date'].split(' - ')[1][:2].replace('/', ''))

            # It would be preferable to use
            # `this_month > last_month`
            # but because of some inconsistency in the data,
            # it is necessary to enforce the stricter constraint
            # of a certain distance away.
            if k > 0 and this_month - last_month > 3:
                year -= 1
            newrow['Year'] = year

            rows.append(newrow)
            k += 1
            
    df = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    title = soup.find('h2', attrs={'class': 'page_title'}).text
    
    return df, title

races = get_races()

# x,y = get_race_data('https://www.realclearpolitics.com/epolls/other/2018_generic_congressional_vote-6185.html')
n = sum([len(r) for r in races.values()])
i = 1
for race_type in races.keys():
    print('{} of {}'.format(i, n))
    for url in races[race_type]:
        df, title = get_race_data(url)
        fn0 = re.sub('[ .]', '_', title)
        fn = 'data/{}-{}.csv'.format(race_type, fn0)
        df.to_csv(fn, index=False)
        i += 1
