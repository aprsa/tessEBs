# Generated by Django 3.0.4 on 2020-09-18 18:40

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0017_auto_20200918_1832'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ephemerissource',
            name='reference',
            field=models.CharField(blank=True, max_length=64, null=True),
        ),
        migrations.AlterField(
            model_name='ephemerissource',
            name='version',
            field=models.CharField(blank=True, max_length=16, null=True),
        ),
    ]
