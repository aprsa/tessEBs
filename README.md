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

```text
# tessEBs/private.py
SECRET_KEY = '[production key as it appears originally in settings.py]'
ENGINE = 'django.db.backends.mysql'  # works for mariadb and mysql
NAME = 'tessEBs'  # name of the mariadb database; has to be tessEBs if settings.py remains unchanged
USER = 'user'  # mariadb username with full access to the database
PASSWORD = 'password'  # mariadb password with full access to the database
ALLOWED_HOSTS = ['allowed_hostname']  # hostname(s) that have access to the database
```

## Running a server:

The current development setup is using django's own server; to run, do (as regular user):
```
screen -S tessEBs
source /path/to/venv/bin/activate
cd /path/to/root/dir
./manage.py runmodwsgi --reload-on-changes
```

## Interfacing the catalog using command line interface (CLI)

The admin interface is fully functional and entries can be added through the admin web form; for bulk ingestion that is of course impractical. To use CLI, use the top-level `manage.py` under the appropriate environment:

```bash
bash$ source /path/to/venv/bin/activate
bash (venv)$ ./manage.py shell
```

### Adding TIC entries to the database

To add a new TIC entry, we can instantiate a `TIC`-class object by passing all relevant parameters:

```python
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

```python
>>> tic = TIC.from_mast(tess_id=1234567890)
>>> tic.save()
```

If a proposed ephemeris is available, the TIC can also be added to triage by using the `add_to_triage()` class method. It takes three arguments: `tess_id`, `ephemeris` and `source`. The `ephemeris` and `source` parameters can be either dictionaries or database model instances. For example:

```python
>>> from catalog.models import TIC, EphemerisSource
>>> tess_id = 7720507
>>> ephemeris = {
...     'bjd0': 1984.533246,
...     'period': 1.577347,
... }
>>> source = EphemerisSource.objects.filter(model='QATS', version='1.0')[0]
>>> tic = TIC.add_to_triage(tess_id=tess_id, ephemeris=ephemeris, source=source)
```

### Creating static files for a binary star

To properly include data files and plots on the website, static files need to be created. To create static files, we use the `create_static_files()` method. The method takes the following arguments:

| Argument | Type | Description |
|----------|------|-------------|
| `static_dir` | string | base directory to store the static files in, 'static/catalog' by default |
| `export_lc` | boolean | should data files be exported to an ascii file |
| `export_spd` | boolean | should a power spectrum be computed and exported to an ascii file |
| `plot_lc` | boolean | should all LC-related plots (lc, zlc, ph) be plotted |
| `plot_spd` | boolean | should a power spectrum plot be plotted |
| `force_overwrite` | boolean | should existing data/plots be overwritten |

Say an EB has already been imported. Static files will be created by issuing:
 
```python
>>> from catalog.models import EB
>>> eb = EB.objects.filter(tic__tess_id=1234567890)[0]
>>> eb.create_static_files(force_overwrite=True)
```
