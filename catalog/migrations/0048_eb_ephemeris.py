# Generated by Django 5.1.6 on 2025-02-26 05:04

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0047_remove_ephemeris_source_ephemeris_source'),
    ]

    operations = [
        migrations.AddField(
            model_name='eb',
            name='ephemeris',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='ephemeris', to='catalog.ephemeris'),
        ),
    ]
