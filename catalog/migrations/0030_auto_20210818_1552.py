# Generated by Django 3.0.4 on 2021-08-18 15:52

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0029_eb_morph_dist'),
    ]

    operations = [
        migrations.AlterField(
            model_name='eb',
            name='date_modified',
            field=models.DateTimeField(auto_now=True, null=True, verbose_name='date modified'),
        ),
    ]
