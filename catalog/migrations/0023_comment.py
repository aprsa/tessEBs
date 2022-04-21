# Generated by Django 3.0.4 on 2020-10-19 18:06

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('catalog', '0022_auto_20201016_1441'),
    ]

    operations = [
        migrations.CreateModel(
            name='Comment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('author', models.CharField(max_length=32, verbose_name='author')),
                ('text', models.TextField(blank=True, null=True, verbose_name='text')),
                ('timestamp', models.DateTimeField(verbose_name='commented on')),
                ('eb', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='comments', to='catalog.EB')),
                ('ephem', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='comments', to='catalog.Ephemeris')),
                ('tic', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='comments', to='catalog.TIC')),
            ],
        ),
    ]
