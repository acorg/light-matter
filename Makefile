.PHONY: check, tcheck, pep8, pyflakes, lint, wc, clean

check:
	python -m discover -v

tcheck:
	trial --rterrors test

pep8:
	find . -name '*.py' -print0 | xargs -0 pep8

pyflakes:
	find . -name '*.py' -print0 | xargs -0 pyflakes

lint: pep8 pyflakes

wc:
	find . -path './light/colors' -prune -o -path './test/colors' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	rm -fr _trial_temp
