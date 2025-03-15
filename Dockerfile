# Use an official Python runtime as a parent image
FROM python:3.12-bookworm

# Install node and npm for bokeh:
RUN apt-get update
RUN apt-get install -y nodejs npm

# Update pip to latest:
RUN pip install --upgrade pip

# Install other packages that the portal needs:
RUN pip install numpy scipy matplotlib astropy astroquery

# Install django stuff:
RUN pip install django django-csv-export-view django-axes mysqlclient uwsgi

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV DJANGO_ALLOWED_HOSTS="localhost 127.0.0.1 [::1]"

# Set the working directory
WORKDIR /opt/build

# Build the custom bokeh version:
COPY deps /opt/build/
RUN pip install bokeh/

# Clean up the build directory:
RUN rm -rf /opt/build/bokeh

# Grab TESS EB server:
RUN git clone https://github.com/aprsa/tessEBs /srv/www/tessEBs

# Open port 8080 to localhost:
EXPOSE 8080
