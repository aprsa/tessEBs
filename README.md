<<<<<<< HEAD
# tessEBs
TESS Eclipsing Binary website
=======
# TESS eclipsing binary catalog (tessEBs)

This code utilizes [django](https://www.djangoproject.com/) to create and maintain the TESS eclipsing binary catalog backbone. The database backend is [mariadb](https://mariadb.org/). It is written and maintained by Andrej Prsa and the TESS Eclipsing Binary Working Group. The associated paper is [ApJS 258, 16, 2022](https://ui.adsabs.harvard.edu/abs/2022ApJS..258...16P/abstract).

### Interfacing the catalog using command line interface (CLI)

The admin interface is fully functional and entries can be added directly; for bulk ingestion that is of course impractical. To use CLI, use the top-level `manage.py` under the appropriate environment:
```
bash$ source /path/to/venv/bin/activate
bash (venv)$ ./manage.py shell
```
### Adding TIC entries to the database
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
Alternatively, we can pull all associated TIC values from MAST:
```
>>> tic = TIC.from_mast(tess_id=1234567890)
>>> tic.save()
```
>>>>>>> 5d1176095aedd2a2d384d52054e0a6712220e882
