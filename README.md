# TESS eclipsing binary catalog (tessEBs)

This code utilizes [django](https://www.djangoproject.com/) to create and maintain the TESS eclipsing binary catalog backbone. The database backend is [mariadb](https://mariadb.org/). It is written and maintained by Andrej Prsa and the TESS Eclipsing Binary Working Group. The associated paper is [ApJS 258, 16, 2022](https://ui.adsabs.harvard.edu/abs/2022ApJS..258...16P/abstract).

## Database structure

The database structure is defined in `catalog/models.py`. The main tables are:

| Table  | Content |
| -----  | ------- |
| TIC    | TESS Input Catalog -- basic properties of the target |
| Sector | TESS sector -- basic sector properties, such as time span and per-camera R.A., Dec and roll |
| Origin | Database entry origin -- who contributed the target candidate to the catalog |
| EB     | Eclipsing Binary entry -- basic properties of the eclipsing binary star |
| Ephemeris | Proposed ephemeris -- eclipse BJD and orbital period |
| EphemerisSource | Source of the proposed ephemeris -- the code used or the author of the manual override |
| Comment | Catalog entry comment -- anything pertinent about the specific eclipsing binary star |

## Installation

The tessEBs backend is developed in python 3.6+ and django 3.0+. The actual deployment versions are python 3.7.6 and django 3.2.13.

To deploy the database (and the website) locally:
* initialize a virtual environment
* install django 3.0+, numpy, mariadb
* create an empty `tessEBs` database
* create `tessEBs/private.py` to set up the server and a link to the database
* run top-level `manage.py` to make all database migrations.

### Setting up credentials

Typically, all settings that pertain to the database and the website reside in `tessEBs/settings.py`. This is certainly a distinct possibility, but given that this file is shared publicly at github, it would be impractical to keep the actual version on github and using a local version in production. Instead, we keep all credentials in the `tessEBs/private.py` file. The minimal contents are:
```
# tessEBs/private.py
SECRET_KEY = '[production key as it appears originally in settings.py]'
ENGINE = 'django.db.backends.mysql'  # works for mariadb and mysql
NAME = 'tessEBs'  # name of the mariadb database; has to be tessEBs if settings.py remains unchanged
USER = 'user'  # mariadb username with full access to the database
PASSWORD = 'password'  # mariadb password with full access to the database
ALLOWED_HOSTS = ['allowed_hostname']  # hostname(s) that have access to the database
```

## Interfacing the catalog using command line interface (CLI)

The admin interface is fully functional and entries can be added through the admin web form; for bulk ingestion that is of course impractical. To use CLI, use the top-level `manage.py` under the appropriate environment:
```
bash$ source /path/to/venv/bin/activate
bash (venv)$ ./manage.py shell
```
### Adding TIC entries to the database

To add a new TIC entry, we can instantiate a `TIC`-class object by passing all relevant parameters:
```
>>> from catalog.models import TIC
>>> tic = TIC(
...     tess_id=1234567890,
...     ra=123.456789,
...     dec=12.3456789,
...     glon=123.456789,
...     glat=12.3456789,
...     Tmag=12.34,
...     teff=12345,
...     logg=3.45,
...     abun=0.123,
...     pmra=0.123,
...     pmdec=0.123,
...     gaia_id=123456789012345
... )
>>> tic.save()
```
Alternatively, we can use a class method to pull all associated TIC values from MAST:
```
>>> tic = TIC.from_mast(tess_id=1234567890)
>>> tic.save()
```
