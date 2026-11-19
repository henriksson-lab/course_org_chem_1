.PHONY: all render preview figures downloads check layout demo-check demo-browser clean

# Keep Quarto's compilation cache inside the writable, ignored project cache.
QUARTO_CACHE_ROOT ?= $(CURDIR)/.quarto/cache

all: render

figures:
	python3 tools/make_pilot_figures.py
	python3 tools/make_foundation_figures.py
	python3 tools/make_reactivity_figures.py
	python3 tools/make_application_figures.py
	python3 tools/make_pacing_foundation_figures.py
	python3 tools/make_pacing_reactivity_figures.py
	python3 tools/make_pacing_application_figures.py
	Rscript tools/make-r-figures.R
	python3 tools/asset_register.py

downloads:
	python3 tools/package_downloads.py

render: downloads
	XDG_CACHE_HOME="$(QUARTO_CACHE_ROOT)" quarto render
	@touch docs/.nojekyll

preview:
	XDG_CACHE_HOME="$(QUARTO_CACHE_ROOT)" quarto preview

check:
	python3 tools/check_course.py
	python3 tools/check_demos.py
	node tools/check_demo_models.js
	python3 tools/check_exercise_code.py
	Rscript instructor/check-calculations.R

layout:
	python3 tools/check_slide_layout.py

demo-check:
	python3 tools/check_demos.py
	node tools/check_demo_models.js

demo-browser:
	python3 tools/check_demo_browser.py

clean:
	rm -rf docs .quarto
