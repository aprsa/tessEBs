# Generated by Django 3.2.13 on 2022-11-05 21:16

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0041_auto_20221104_2306'),
    ]

    operations = [
        migrations.AlterField(
            model_name='tic',
            name='asas_id',
            field=models.CharField(blank=True, max_length=20, null=True),
        ),
    ]
