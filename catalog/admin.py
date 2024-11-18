from django.contrib import admin

from .models import TIC, EB, Origin, Sector, EphemerisSource, Ephemeris, Comment, BLS_run

#class EBAdmin(admin.ModelAdmin):
#    model = EB
#    filter_horizontal = ('source',)


class EphemerisInline(admin.StackedInline):
    model = Ephemeris
    extra = 0


class TICAdmin(admin.ModelAdmin):
    model = TIC
    inlines = [EphemerisInline]


class EBAdmin(admin.ModelAdmin):
    model = EB
    # inlines = [EphemerisInline]
    search_fields = ('tic__tess_id',)
    # list_display = ('ephemerides',)


admin.site.register(TIC, TICAdmin)
admin.site.register(EB, EBAdmin)
admin.site.register(Sector)
admin.site.register(Origin)
admin.site.register(EphemerisSource)
admin.site.register(Ephemeris)
admin.site.register(Comment)
admin.site.register(BLS_run)
