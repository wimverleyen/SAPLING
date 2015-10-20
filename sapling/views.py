from __future__ import absolute_import

from django.core.urlresolvers import reverse_lazy
from django.http import HttpResponse, HttpResponseRedirect
from django.contrib import messages
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.views.generic import TemplateView, CreateView, FormView, RedirectView
from sapling.settings import BASE_DIR

from braces.views import AnonymousRequiredMixin, LoginRequiredMixin, FormValidMessageMixin, MessageMixin

from .forms import RegistrationForm, LoginForm


def send_file1(request):
  import os
  from django.core.servers.basehttp import FileWrapper
  import mimetypes

  filename = BASE_DIR+"/sapling/media/sapling_sample_attention-deficit-hyperactivity-disorder_17.pdf"
  download_name = "sapling_sample_attention-deficit-hyperactivity-disorder_17.pdf"
  wrapper = FileWrapper(open(filename))
  content_type = mimetypes.guess_type(filename)[0]
  response = HttpResponse(wrapper,content_type=content_type)
  response['Content-Length'] = os.path.getsize(filename)
  response['Content-Disposition'] = "attachment; filename=%s"%download_name
  return response


def send_file2(request):
  import os
  from django.core.servers.basehttp import FileWrapper
  import mimetypes

  filename = BASE_DIR+"/sapling/media/sapling_sample_iossifov_asd_probands_recurrent_2014_1.pdf"
  download_name = "sapling_sample_iossifov_asd_probands_recurrent_2014_1.pdf"
  wrapper = FileWrapper(open(filename))
  content_type = mimetypes.guess_type(filename)[0]
  response = HttpResponse(wrapper,content_type=content_type)
  response['Content-Length'] = os.path.getsize(filename)
  response['Content-Disposition'] = "attachment; filename=%s"%download_name
  return response


def send_file3(request):
  import os
  from django.core.servers.basehttp import FileWrapper
  import mimetypes

  filename = BASE_DIR+"/sapling/media/sapling_sample_synsysnet_4_22_2013_8.pdf"
  download_name = "sapling_sample_synsysnet_4_22_2013_8.pdf"
  wrapper = FileWrapper(open(filename))
  content_type = mimetypes.guess_type(filename)[0]
  response = HttpResponse(wrapper,content_type=content_type)
  response['Content-Length'] = os.path.getsize(filename)
  response['Content-Disposition'] = "attachment; filename=%s"%download_name
  return response


class HomePageView(TemplateView):
  template_name = "home.html"


class AboutPageView(TemplateView):
  template_name = "about.html"


class ContactPageView(TemplateView):
  template_name = "contact.html"


class SignUpView(AnonymousRequiredMixin, FormValidMessageMixin, CreateView):
  form_class = RegistrationForm
  form_valid_message = "Thanks for signing up! Go ahead and login"
  model = User
  success_url = reverse_lazy('login')
  template_name = "accounts/signup.html"


class LoginView(AnonymousRequiredMixin, FormValidMessageMixin, FormView):
  form_class = LoginForm
  form_valid_message = "You are logged in now!"
  success_url = reverse_lazy('home')
  template_name = "accounts/login.html"

  def form_valid(self, form):
    username = form.cleaned_data['username']
    password = form.cleaned_data['password']

    user = authenticate(username=username, password=password)

    if user is not None and user.is_active:
      login(self.request, user)
      return super(LoginView, self).form_valid(form)
    else:
      return self.form_invalid(form)


class LogoutView(LoginRequiredMixin, MessageMixin, RedirectView):
  url = reverse_lazy('home')

  def get(self, request, *args, **kwargs):
    messages.success(request, "You have been logged out")
    print "logout user view called"
    logout(request)
    return super(LogoutView, self).get(request, *args, **kwargs)


def custom_404(request):
  return render_to_response('404.html')


def custom_500(request):
  return render_to_response('500.html')
