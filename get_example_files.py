import urllib
import os.path

if not os.path.isfile('flash.p'):
    urllib.urlretrieve ('https://s3.amazonaws.com/pybloch/flash.p', 'flash.p')

if not os.path.isfile('ssfp_ss_startup.p'):
    urllib.urlretrieve ('https://s3.amazonaws.com/pybloch/ssfp_ss_startup.p', 'ssfp_ss_startup.p')

if not os.path.isfile('ssfp_store.p'):
    urllib.urlretrieve ('https://s3.amazonaws.com/pybloch/ssfp_store.p', 'ssfp_store.p')
