# Generated by Django 3.0.4 on 2020-07-23 14:54

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0014_auto_20200722_1548'),
    ]

    operations = [
        migrations.AlterField(
            model_name='eb',
            name='sectors',
            field=models.ManyToManyField(blank=True, null=True, to='catalog.Sector'),
        ),
    ]
