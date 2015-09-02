from django.conf import settings
from django.conf.urls.static import static
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf.urls import patterns, include, url

from .views import HomePageView, SignUpView, LoginView, LogoutView, AboutPageView, ContactPageView, send_file1, send_file2, send_file3


handler_404 = 'saplingdev.views.custom_404'
handler_500 = 'saplingdev.views.custom_500'

urlpatterns = patterns('',
  url(regex = r'^$',
      view = HomePageView.as_view(),
      name = 'home'),

  url(r'^media/sapling_sample_attention-deficit-hyperactivity-disorder_17.pdf', send_file1),
  url(r'^media/sapling_sample_iossifov_asd_probands_recurrent_2014_1.pdf', send_file2),
  url(r'^media/sapling_sample_synsysnet_4_22_2013_8.pdf', send_file3),

  url(r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root':settings.STATIC_ROOT}),

  url(regex=r'accounts/register/$', 
      view=SignUpView.as_view(),
      name='signup'),
  url(regex=r'accounts/login/$', 
      view=LoginView.as_view(),
      name='login'),
  url(regex=r'accounts/logout/$', 
      view=LogoutView.as_view(),
      name='logout'),

  url(regex=r'^resetpassword/$', 
      view='django.contrib.auth.views.password_reset', 
      name='password_reset'),
  url(regex=r'^resetpassword/passwordsent/$', 
      view='django.contrib.auth.views.password_reset_done', 
      name='password_reset_done'),
  url(regex=r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>.+)/$', 
      view='django.contrib.auth.views.password_reset_confirm', 
      name='password_reset_confirm'),
  url(regex=r'^reset/done/$', 
      view='django.contrib.auth.views.password_reset_complete', 
      name='password_reset_complete'),

  url(regex = r'^about/$',
      view = AboutPageView.as_view(),
      name = 'about'),
  url(regex = r'^contact/$',
      view = ContactPageView.as_view(),
      name = 'contact'),

  url(r'^gfps/', include('GFPdevs.urls', namespace='gfps'))
)

if settings.DEBUG:
  urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
