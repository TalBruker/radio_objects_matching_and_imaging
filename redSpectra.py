import numpy as np
from astropy import units as u
import astropy.cosmology as cosmo
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from PIL import Image
import requests
from io import BytesIO
import math
from astropy.io import fits
from scipy.optimize import curve_fit


def excelToPictures(excelName, starName, raName, decName, saveName, filters="grizy", redName=None):
    data = pd.read_csv(r'{}'.format(excelName))
    print(data)
    names = pd.DataFrame(data, columns=['{}'.format(starName)])
    print(names)
    nameList = data['{}'.format(starName)].drop_duplicates().to_list()
    print(nameList)
    RaList = data['{}'.format(raName)].drop_duplicates().to_list()
    DecList = data['{}'.format(decName)].drop_duplicates().to_list()
    if redName is not None:
        RedList = data['{}'.format(redName)].to_list()
    else:
        RedList = np.empty(len(RaList))
        RedList[:] = np.NaN
    print('The redshift array is: ', RedList)

    def getimages(ra, dec, name=None, size=240, filters="grizy"):
        """Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
               "&filters={filters}").format(**locals())
        if (name is not None):
            url = ("{service}?skycekk={name}&size={size}&format=fits"
                   "&filters={filters}").format(**locals())
        table = Table.read(url, format='ascii')
        return table

    def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
        """Get URL for images in the table

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png" or "fits")
        color = if True, creates a color image (only for jpg or png format).
                Default is return a list of URLs for single-filter grayscale images.
        Returns a string with the URL
        """

        if color and format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = getimages(ra, dec, size=size, filters=filters)
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
        if output_size:
            url = url + "&output_size={}".format(output_size)
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[np.argsort(flist)]
        if color:
            if len(table) > 3:
                # pick 3 filters
                table = table[[0, len(table) // 2, len(table) - 1]]
            for i, param in enumerate(["red", "green", "blue"]):
                url = url + "&{}={}".format(param, table['filename'][i])
        else:
            urlbase = url + "&red="
            url = []
            for filename in table['filename']:
                url.append(urlbase + filename)
        return url

    def getcolorim(ra, dec, size=240, output_size=None, filters="grizy", format="jpg"):
        """Get color image at a sky position

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png")
        Returns the image
        """

        if format not in ("jpg", "png"):
            raise ValueError("format must be jpg or png")
        url = geturl(ra, dec, size=size, filters=filters, output_size=output_size, format=format, color=True)
        r = requests.get(url)
        im = Image.open(BytesIO(r.content))
        return im

    def getgrayim(ra, dec, size=240, output_size=None, filter="g", format="jpg"):
        """Get grayscale image at a sky position

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filter = string with filter to extract (one of grizy)
        format = data format (options are "jpg", "png")
        Returns the image
        """

        if format not in ("jpg", "png"):
            raise ValueError("format must be jpg or png")
        if filter not in list("grizy"):
            raise ValueError("filter must be one of grizy")
        url = geturl(ra, dec, size=size, filters=filter, output_size=output_size, format=format)
        r = requests.get(url[0])
        im = Image.open(BytesIO(r.content))
        return im

    ra = RaList[0]
    dec = DecList[0]
    imsize = 60 * u.arcsecond
    ImList = getcolorim(ra, dec, filters=filters)
    if (ImList is None):
        print('the image for the galaxi {} is out of boundry'.format(starName[0]))
    plt.imshow(ImList)
    size = (imsize * cosmo.WMAP9.kpc_proper_per_arcmin(RedList[0])).decompose().to(u.kpc)
    plt.title(f"{size:0.03f}")
    linestart = (10, 10)  # the size of the picture is 235 pixels.
    pixelsize = 235 * 10 * u.kpc
    lineEnd = (10, 10 + pixelsize.value / size.value)  # (10,10 + size.value/2)
    plt.plot(linestart, lineEnd, linewidth=10, color='r')
    plt.show()
    pylab.rcParams.update({'font.size': 12})

    print("size: ", len(RaList))

    size = (imsize * cosmo.WMAP9.kpc_proper_per_arcmin(0.130635)).decompose().to(u.kpc)
    print(f"{size:0.03f}")

    fig = plt.figure()

    outOfRangeStars = []

    columns = 10
    rows = 10
    loop_time = int((len(DecList) - 1) / 100) + 1
    for i in range(loop_time):
        fig = plt.figure()
        for j in range(columns * rows):  # columns*rows):
            index = 100 * i + j
            try:
                img = getcolorim(RaList[index], DecList[index], filters=filters)
            except IndexError:
                outOfRangeStars.append(nameList)
                continue
            fig.add_subplot(rows, columns, j + 1)
            plt.imshow(img)
            plt.axis('off')
            size = (imsize * cosmo.WMAP9.kpc_proper_per_arcmin(RedList[index])).decompose().to(u.kpc)
            linestart = (24, 24)
            lineEnd = (24, 24 + pixelsize.value / size.value)
            plt.plot(linestart, lineEnd, linewidth=1, color='r')
            # plt.savefig('galaxies.svg', format='svg')
        plt.savefig('{}_{}.png'.format(saveName, i), dpi=1000)
        plt.show()


#    rows = int(0.9 + len(DecList) / 10)
#    for i in range(len(DecList)):
#        try:
#            img = getcolorim(RaList[i], DecList[i], filters=filters)
#        except IndexError:
#            outOfRangeStars.append(nameList)
#            continue
#        fig.add_subplot(rows, columns, i + 1)
#        plt.imshow(img)
#        plt.axis('off')
#        size = (imsize * cosmo.WMAP9.kpc_proper_per_arcmin(RedList[i])).decompose().to(u.kpc)
#        linestart = (24, 24)
#        lineEnd = (24, 24 + pixelsize.value / size.value)
#        plt.plot(linestart, lineEnd, linewidth=1, color='r')

# plt.savefig('galaxies.svg', format='svg')
# print('The stars out of the angular range are: ', outOfRangeStars)


# excelToPictures('Astro-project.csv','Name','RA','DEC','galaxies1color','Redshift')
# excelToPictures('blazars.csv','Identifier','RA (ICRS (J2000))','DEC (ICRS (J2000))','blazarGalaxiesColorIm')
# excelToPictures('blazars.csv','Identifier','RA (ICRS (J2000))','DEC (ICRS (J2000))','blazarGalaxiesGZYIm',filters="yzg")
# excelToPictures('edgeOnGalaxiesCheck.csv','specobjid','ra','dec','edgeOnPicCheck',redName='redshift')
# excelToPictures('faceOnGalaxiesCheck.csv','specobjid','ra','dec','faceOnPicCheck',redName='redshift')

from astropy.table import Table
from astropy.coordinates import SkyCoord

import csv


def ra_to_csv(red_data, name):
    #    data_angle = red_data[['ra','dec']]
    #    print(data_angle)
    #    data_angle.to_csv('red_edge_on_galaxies_placemnent.csv')
    with open(name, 'w', newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter=',')
        my_writer.writerows(red_data)


def data_matching():
    dat = Table.read('first_14dec17.fits', format='fits')
    df = dat.to_pandas()
    radio_obj_ra = np.array(df['RA'].to_list())
    radio_obj_dec = np.array(df['DEC'].to_list())
    radio_coord = SkyCoord(ra=radio_obj_ra * u.degree, dec=radio_obj_dec * u.degree)
    data = pd.read_csv('edge_On_Galaxies_With_Average_Redshift_talbruker.csv')
    data_ra = np.array(data['ra'].to_list())
    data_dec = np.array(data['dec'].to_list())
    data_coord = SkyCoord(ra=data_ra * u.degree, dec=data_dec * u.degree)
    max_sep = 1.0 * u.arcmin
    idx, d2d, d3d = data_coord.match_to_catalog_sky(radio_coord)

    def plot_number_of_stars():
        def number_of_matchs(max_sep):
            sep_constraint = d2d < max_sep
            data_matches = data_coord[sep_constraint]
            # radio_coord_matches = radio_coord[idx[sep_constraint]]
            return (len(data_matches))

        X = np.linspace(0, 1, 10000)

        Y = [number_of_matchs(x * u.arcmin) for x in zip(X.tolist())]
        plt.plot(X, Y)
        plt.show()

        X = np.linspace(0, 0.1, 10000)

        Y = [number_of_matchs(x * u.arcmin) for x in zip(X.tolist())]
        plt.plot(X, Y)
        plt.show()

        def parab(x, a, c):
            return a + c * x ** 2  # + b*x

        X = np.linspace(0.1, 1, 1000)
        Y = [number_of_matchs(x * u.arcmin) for x in zip(X.tolist())]
        popt, pcov = curve_fit(parab, X, Y)
        print('popt is : _____________________________________')
        print(popt)
        X = np.linspace(0, 1, 10000)
        Y = [number_of_matchs(x * u.arcmin) for x in zip(X.tolist())]
        plt.plot(X, Y)
        plt.plot(X, parab(X, *popt))
        plt.xlabel(r'$\theta$ [arcminutes]')
        plt.ylabel('Number Of Matching Stars')
        plt.title(r'The Number Of Matching Galaxies As a Function Of The' +'\n' +r'Separation Angle $\theta$')
        plt.savefig('noise_in_galaxy_num.pdf')
        plt.show()
        print(popt)

    #plot_number_of_stars()


    angle_sep = 0.04 * u.arcmin
    sep_constraint = d2d < angle_sep
    data_matches = data_coord[sep_constraint]
    catalog_matches = radio_coord[idx[sep_constraint]]

    data_matches_array = [i.split() for i in data_matches.to_string()]

    return data_matches


def fit_names(data_matches, secondary = None):
    def string_cut(string, cut_string, before, until=0, from_start=False):
        if from_start:
            return string[0:string.find(cut_string) + until]
        else:
            return string[string.find(cut_string) + len(cut_string) - 1 - before:string.find(cut_string) + until]

    #print('__________________________________________')
    red_data_ra = []
    red_data_dec = []
    # checking the cose for first 20 galaxies
    for i in data_matches:
        rahmsstr = i.ra.to_string(u.hour)
        decdmsstr = i.dec.to_string(u.degree, alwayssign=True)
        #print(rahmsstr, decdmsstr)
        hour_string = string_cut(rahmsstr, 'h', 0, 0, True)
        if len(hour_string) == 1:
            hour_string = '0' + hour_string
        minute_string = string_cut(rahmsstr, 'm', 2, 0)
        minute_string = minute_string + '0'
        red_data_ra.append(int(hour_string + minute_string))

        d_string = string_cut(decdmsstr, 'd', 0, 0, True)
        d_minute_string = string_cut(decdmsstr, 'm', 2, 0)
        d_second_string = string_cut(decdmsstr, 'm', -1, 3)
        d_minute_string = d_minute_string + '0'
        red_data_dec.append(int(d_string + d_minute_string))

    #print(red_data_ra)
    #print(red_data_dec)

    def print_fits():
        from astropy.io import fits
        hdu_list = fits.open('coverage-north-3arcmin-14dec17.fits')
        hdu_list.info()
        print(type(hdu_list))
        image_data = fits.getdata('firstcutout _1__1625665026')

        plt.figure()
        plt.imshow(image_data, cmap='gray')
        plt.colorbar()
        plt.show()

    url = 'https://archive.stsci.edu/pub/vla_first/data/'

    # rahmsstr[rahmsstr.find(' ') + len(' ') + 1:rahmsstr.find(' ') + len(' ') + 4]

    def listFD(url):
        page = requests.get(url).text
        page = page.rstrip().split('\n')
        page = [i[i.find('<a href="') + len('<a href="'):i.find('/">')] for i in page]
        ra = []
        for i in page:
            if i.isdecimal():
                ra.append(int(i))
        return ra

    def find_nearest(array, value, secondary = None):
        darray = np.asarray(array)
        idx = (np.abs(darray - value)).argmin()
        maxim = (np.abs(darray - value)).argmax()
        if secondary is not None:
            for i in range(secondary):
                array[idx] = array[maxim]
                darray = np.asarray(array)
                idx = (np.abs(darray - value)).argmin()
        return darray[idx]

    ra_galaxies_array = [find_nearest(listFD(url), i,secondary) for i in red_data_ra]
    # now we have the ra name of the galaxy fit

    for i in range(len(ra_galaxies_array)):
        ra_galaxies_array[i] = str(ra_galaxies_array[i])
        zero_add = 5 - len(ra_galaxies_array[i])
        for j in range(zero_add):
            ra_galaxies_array[i] = '0' + ra_galaxies_array[i]

    def listFD2(url):
        page = requests.get(url).text
        page = page.rstrip().split('\n')
        page = [i[i.find('.fits">') + len('.fits">') + 5:i.find('.fits</a>') - 1] for i in page]
        ra = []
        ra_index = -1
        for i in page:
            if len(i) > 0:
                if i[0] == '+':
                    if i[1:].isdecimal():
                        if ra_index == -1 or int(i[1:]) != ra[ra_index]:
                            ra.append(int(i[1:]))
                            ra_index += 1
                if i[0] == '-':
                    if i[1:].isdecimal():
                        if ra_index == -1 or int(i[1:]) != -ra[ra_index]:
                            ra.append(-int(i[1:]))
                            ra_index += 1
        return ra

    #print('_____________________________________')
    dec_galaxies_array = []
    for i in range(len(ra_galaxies_array)):
        url_new = url + ra_galaxies_array[i] + '/'
        dec_galaxies_array.append(find_nearest(listFD2(url_new), red_data_dec[i]))

    #print(dec_galaxies_array)

    for i in range(len(dec_galaxies_array)):
        string = str(dec_galaxies_array[i])
        if string[0] == '-':
            zero_add = 5 - (len(string) - 1)
            if zero_add > 0:
                string = list(string)
                string[0] = '0'
                string = ''.join(string)
                for j in range(zero_add - 1):
                    string = '0' + string
                string = '-' + string
        else:
            zero_add = 5 - len(string)
            for j in range(zero_add):
                string = '0' + string
            string = '+' + string
        dec_galaxies_array[i] = string

    name_array = []

    for i in range(len(dec_galaxies_array)):
        name_array.append(ra_galaxies_array[i] + dec_galaxies_array[i])
    name_array = [[i] for i in name_array]

    #ra_to_csv(name_array, 'red_edge_on_galaxies_fits_name.csv')
    return name_array


location_list = pd.read_csv(r'red_edge_on_galaxies_fits_name.csv').values.tolist()
location_list = [i[0] for i in location_list]

a = data_matching()

import requests
from shutil import copyfile
import urllib.request
from astropy.utils.data import download_file, clear_download_cache
from astropy.config.paths import get_cache_dir
from astropy.wcs import WCS


# https://forum.qiime2.org/t/no-space-left-on-device-for-qiime-feature-classifier/2122
def get_fit(fit_name, star_loc):
    ra_dec = [star_loc.ra.degree, star_loc.dec.degree]
    star_location = [SkyCoord(ra_dec[0] * u.degree - 1 * u.arcmin, ra_dec[1] * u.degree - 1 * u.arcmin),
                   SkyCoord(ra_dec[0] * u.degree + 1 * u.arcmin, ra_dec[1] * u.degree + 1 * u.arcmin)]
    url = 'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5]
    page = requests.get(url).text
    page = page.rstrip().split('\n')
    page = [i[i.find('.fits">') + len('.fits">'):i.find('.fits</a>')] for i in page]
    name_array = [i for i in page if fit_name in i]
    fit_file = name_array[len(name_array) - 1]
    fit_url = r'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5] + r'/' + fit_file + r'.fits'
    # r = requests.get(fit_url, allow_redirects=True)
    # src = url
    # dst = r'C:\Users\tal2\Desktop\web scraping\temp'
    #    with open(r'C:\Users\tal2\Desktop\web scraping\temp\temp_fit.fit', 'wb') as f:
    #        f.write(r.content)

    image_file = download_file(fit_url, cache=True)
    hdu_list = fits.open(image_file)
    # hdu_list.info()
    header = hdu_list[0].header
    image_data = hdu_list[0].data
    image_data = image_data[0, 0, :, :]
    header['NAXIS'] = 2
    for i in range(3, 5):
        for j in ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']:
            string = j + str(i)
            header.remove(string, ignore_missing=True, remove_all=True)
    hdu_list.close()
    wcs_helix = WCS(header)
    # (200.06152002, 2.6541609)
    a = wcs_helix.world_to_pixel(star_location[0])
    b = wcs_helix.world_to_pixel(star_location[1])

    def np_tup_to_int(tup):
        return [round(tup[1].tolist()), round(tup[0].tolist())]

    a = np_tup_to_int(a)
    b = np_tup_to_int(b)
    c = [min(a[0], b[0]), min(a[1], b[1])]
    d = [max(a[0], b[0]), max(a[1], b[1])]

    size = image_data.shape
    if d[0] > size[0] or d[1] > size[1]:
        fit_name = fit_names([star_loc], secondary=True)[0][0]
        url = 'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5]
        page = requests.get(url).text
        page = page.rstrip().split('\n')
        page = [i[i.find('.fits">') + len('.fits">'):i.find('.fits</a>')] for i in page]
        name_array = [i for i in page if fit_name in i]
        fit_file = name_array[len(name_array) - 1]
        fit_url = r'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5] + r'/' + fit_file + r'.fits'
        image_file = download_file(fit_url, cache=True)
        hdu_list = fits.open(image_file)
        # hdu_list.info()
        header = hdu_list[0].header
        image_data = hdu_list[0].data
        image_data = image_data[0, 0, :, :]
        header['NAXIS'] = 2
        for i in range(3, 5):
            for j in ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']:
                string = j + str(i)
                header.remove(string, ignore_missing=True, remove_all=True)
        hdu_list.close()
        wcs_helix = WCS(header)
        # (200.06152002, 2.6541609)
        a = wcs_helix.world_to_pixel(star_location[0])
        b = wcs_helix.world_to_pixel(star_location[1])

        a = np_tup_to_int(a)
        b = np_tup_to_int(b)
        c = [min(a[0], b[0]), min(a[1], b[1])]
        d = [max(a[0], b[0]), max(a[1], b[1])]

    return [image_data[c[0]:d[0], c[1]:d[1]], url]
    def comment():
        if d[0] > size[0] or d[1] > size[1]:
                fit_name = fit_names([star_loc],secondary = True)[0][0]
                url = 'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5]
                page = requests.get(url).text
                page = page.rstrip().split('\n')
                page = [i[i.find('.fits">') + len('.fits">'):i.find('.fits</a>')] for i in page]
                name_array = [i for i in page if fit_name in i]
                fit_file = name_array[len(name_array) - 1]
                fit_url = r'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5] + r'/' + fit_file + r'.fits'
                image_file = download_file(fit_url, cache=True)
                hdu_list = fits.open(image_file)
                # hdu_list.info()
                header = hdu_list[0].header
                image_data = hdu_list[0].data
                image_data = image_data[0, 0, :, :]
                header['NAXIS'] = 2
                for i in range(3, 5):
                    for j in ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']:
                        string = j + str(i)
                        header.remove(string, ignore_missing=True, remove_all=True)
                hdu_list.close()
                wcs_helix = WCS(header)
                # (200.06152002, 2.6541609)
                a = wcs_helix.world_to_pixel(star_location[0])
                b = wcs_helix.world_to_pixel(star_location[1])


                a = np_tup_to_int(a)
                b = np_tup_to_int(b)
                c = [min(a[0], b[0]), min(a[1], b[1])]
                d = [max(a[0], b[0]), max(a[1], b[1])]



    # fig = plt.figure(figsize=(10, 10))
    # ax = plt.subplot(projection=wcs_helix)
    # plt.imshow(image_data[c[0]:d[0],c[1]:d[1]], cmap = 'gray')
    # plt.xlabel(r'RA')
    # plt.ylabel(r'Dec')
    # plt.show()
    return [image_data[c[0]:d[0], c[1]:d[1]],fit_name]

def get_fit_from_full_name(fit_name, star_loc):
    ra_dec = [star_loc.ra.degree, star_loc.dec.degree]
    star_location = [SkyCoord(ra_dec[0] * u.degree - 1 * u.arcmin, ra_dec[1] * u.degree - 1 * u.arcmin),
                   SkyCoord(ra_dec[0] * u.degree + 1 * u.arcmin, ra_dec[1] * u.degree + 1 * u.arcmin)]

    fit_url = r'https://archive.stsci.edu/pub/vla_first/data/' + fit_name[0:5] + r'/' + fit_name + r'.fits'
    # r = requests.get(fit_url, allow_redirects=True)
    # src = url
    # dst = r'C:\Users\tal2\Desktop\web scraping\temp'
    #    with open(r'C:\Users\tal2\Desktop\web scraping\temp\temp_fit.fit', 'wb') as f:
    #        f.write(r.content)

    image_file = download_file(fit_url, cache=True)
    hdu_list = fits.open(image_file)
    # hdu_list.info()
    header = hdu_list[0].header
    image_data = hdu_list[0].data
    image_data = image_data[0, 0, :, :]
    header['NAXIS'] = 2
    for i in range(3, 5):
        for j in ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']:
            string = j + str(i)
            header.remove(string, ignore_missing=True, remove_all=True)
    hdu_list.close()
    wcs_helix = WCS(header)
    # (200.06152002, 2.6541609)
    a = wcs_helix.world_to_pixel(star_location[0])
    b = wcs_helix.world_to_pixel(star_location[1])

    def np_tup_to_int(tup):
        return [round(tup[1].tolist()), round(tup[0].tolist())]

    a = np_tup_to_int(a)
    b = np_tup_to_int(b)
    c = [min(a[0], b[0]), min(a[1], b[1])]
    d = [max(a[0], b[0]), max(a[1], b[1])]

    # fig = plt.figure(figsize=(10, 10))
    # ax = plt.subplot(projection=wcs_helix)
    # plt.imshow(image_data[c[0]:d[0],c[1]:d[1]], cmap = 'gray')
    # plt.xlabel(r'RA')
    # plt.ylabel(r'Dec')
    # plt.show()
    return [image_data[c[0]:d[0], c[1]:d[1]],fit_name]

#excelToPictures('red_edge_on_galaxies_placemnent_code.csv', 'RA', 'RA', 'DEC', 'galaxies_edge_with_radio_emission')

def print_fit(fits_name_list):
    matched_data = data_matching()
    print('len matched data is:', len(matched_data))
    arr = []
    loop_time = int((len(matched_data) - 1) / 100) + 1

    for i in [9]:#range(loop_time):
        print(i)
        columns = 10
        rows = 10
        fig = plt.figure()
        for j in range(columns*rows):
            index = 100 * i + j
            fits_file_name = fits_name_list[index]
            star_loc = matched_data[index]
            #ra_dec = [matched_data[index].ra.degree, matched_data[index].dec.degree]
            #data_coords = [SkyCoord(ra_dec[0] * u.degree - 1 * u.arcmin, ra_dec[1] * u.degree - 1 * u.arcmin),
            #               SkyCoord(ra_dec[0] * u.degree + 1 * u.arcmin, ra_dec[1] * u.degree + 1 * u.arcmin)]
            fig.add_subplot(rows, columns, j + 1)
            img, new = get_fit_from_full_name(fits_file_name, star_loc)
            arr.append([new])
            plt.imshow(np.fliplr(img))
            plt.axis('off')
            #try:
            #    plt.contour(np.fliplr(img), levels = 3,colors='white', alpha=0.3)
            #except:
            #    print("error occured in galaxy number {}".format(index + 1))

            if index == len(fits_name_list)-1:
                break


        clear_download_cache()

        # plt.savefig('galaxies.svg', format='svg')
        plt.savefig('radio_objects_{}.png'.format(i), dpi=1000)
        plt.show()

    def comment():
        for i in range(len([254, 258, 536])):
            index = [254, 258, 536][i]
            columns = 1
            rows = 3
            fig = plt.figure()
            fits_file_name = fits_name_list[i]
            star_loc = matched_data[index]
            # ra_dec = [matched_data[index].ra.degree, matched_data[index].dec.degree]
            # data_coords = [SkyCoord(ra_dec[0] * u.degree - 1 * u.arcmin, ra_dec[1] * u.degree - 1 * u.arcmin),
            #               SkyCoord(ra_dec[0] * u.degree + 1 * u.arcmin, ra_dec[1] * u.degree + 1 * u.arcmin)]
            fig.add_subplot(rows, columns, i + 1)
            img, new = get_fit(fits_file_name, star_loc)
            arr.append([new])
            plt.imshow(np.fliplr(img))
            plt.axis('off')
            plt.show()


    #ra_to_csv(arr, 'red_edge_on_galaxies_fits_name_new.csv')



def catalog_information(characteristic):
    #characteristic is a string from the list: 'RA', 'DEC', 'SIDEPROB', 'FPEAK', 'FINT', 'RMS', 'MAJOR', 'MINOR',
    #   'POSANG', 'FITTED_MAJOR', 'FITTED_MINOR', 'FITTED_POSANG', 'FLDNAME',
    #   'NSDSS', 'SDSS_SEP', 'SDSS_MAG', 'SDSS_CLASS', 'NTMASS', 'TMASS_SEP',
    #   'TMASS_MAG', 'YEAR', 'MJD', 'MJDRMS', 'MJDSTART', 'MJDSTOP'
    dat = Table.read('first_14dec17.fits', format='fits')
    df = dat.to_pandas()
    radio_obj_ra = np.array(df['RA'].to_list())
    radio_obj_dec = np.array(df['DEC'].to_list())
    radio_coord = SkyCoord(ra=radio_obj_ra * u.degree, dec=radio_obj_dec * u.degree)
    data = pd.read_csv('edge_On_Galaxies_With_Average_Redshift_talbruker.csv')
    data_ra = np.array(data['ra'].to_list())
    data_dec = np.array(data['dec'].to_list())
    data_coord = SkyCoord(ra=data_ra * u.degree, dec=data_dec * u.degree)
    idx, d2d, d3d = data_coord.match_to_catalog_sky(radio_coord)


    angle_sep = 0.04 * u.arcmin
    sep_constraint = d2d < angle_sep
    data_matches = data_coord[sep_constraint]
    catalog_matches = radio_coord[idx[sep_constraint]]


    print(data_matches)
    print(catalog_matches)
    return [np.array(df[i].to_list())[idx[sep_constraint]] for i in characteristic]


#print(catalog_information(['FPEAK','FINT']))



data = catalog_information(['FINT','SDSS_MAG'])
f_hiluk = np.log10(data[0])
for i in range(len(data[0])):
    f_hiluk[i] = f_hiluk[i] - (data[1][i] - 8.9 + 6)/2.5

good_stars = []
bad_stars = []
for i in range(len(f_hiluk)):
    if f_hiluk[i] > -15:
        good_stars.append(f_hiluk[i])
    else:
        bad_stars.append(i)

#print(bad_stars) [254, 258, 536]
loc_bad_stars = []
for i in bad_stars:
    loc_bad_stars.append(location_list[i])

fit_name_bad_stars = catalog_information(['FLDNAME'])[0]
fit_name_bad_stars = [i.decode('UTF-8')[0:-1] for i in fit_name_bad_stars]
#print(fit_name_bad_stars)
fit_name_bad_stars = np.array(fit_name_bad_stars)[bad_stars]

#print_fit(fit_name_bad_stars)

plt.hist(good_stars,bins = 40)

plt.xlabel(r'$log_{10}\left(\frac{F_{\nu,radio}}{F_{\nu,optical}}\right)$')

plt.title(r'histogram of the galaxies as function of their' +'\n' + r'$log_{10}\left(\frac{F_{\nu,radio}}{F_{\nu,optical}}\right)$')
#plt.savefig('histogram_of_stars.png', dpi=1000)
plt.show()

# No of data points used
f_hiluk = np.array(f_hiluk)
N = len(f_hiluk)

# sort the data in ascending order
x = np.sort(f_hiluk)

cut_num = 10
N = N- cut_num
x = x[cut_num:]

# normal distribution
data = np.random.randn(N)



# get the cdf values of y
y = np.arange(N) / float(N)

# plotting
plt.xlabel(r'$log_{10}\left(\frac{F_{\nu,radio}}{F_{\nu,optical}}\right)$')
plt.ylabel('cumulative distribution function')

plt.title(r'cumulative distribution function of the galaxies as function of their' +'\n' + r'$log_{10}\left(\frac{F_{\nu,radio}}{F_{\nu,optical}}\right)$')

plt.plot(x, y)
#plt.savefig('cdf_of_stars.png', dpi=1000)
plt.show()

#print(location_list[14])
#ra_dec = [a[14].ra.degree, a[14].dec.degree]
#data_coords = [SkyCoord(ra_dec[0] * u.degree - 1 * u.arcmin, ra_dec[1] * u.degree - 1 * u.arcmin),
#               SkyCoord(ra_dec[0] * u.degree + 1 * u.arcmin, ra_dec[1] * u.degree + 1 * u.arcmin)]

#img = get_fit('17210+58155', data_coords)
#plt.imshow(img)
#plt.show()

#print('_________________')
#print(fit_names([a[14]],secondary = True)[0][0])
#print_fit()
#data_matching()
fit_name = catalog_information(['FLDNAME'])[0]
fit_name = [i.decode('UTF-8') for i in fit_name]
print(fit_name)
print(len(fit_name))

#print(location_list[0])
print_fit(fit_name)
print(len(f_hiluk))