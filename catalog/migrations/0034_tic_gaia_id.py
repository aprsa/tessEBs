# Generated by Django 3.0.4 on 2021-10-28 00:25

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0033_auto_20210916_1442'),
    ]

    operations = [
        migrations.AddField(
            model_name='tic',
            name='gaia_id',
            field=models.BigIntegerField(blank=True, null=True, verbose_name='gaia id'),
        ),
    ]