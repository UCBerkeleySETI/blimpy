### Contributing to blimpy

Thanks for thinking about contributing code to `blimpy`! Here's a few quick guidelines:

* You'll need to follow our [code of conduct](https://github.com/UCBerkeleySETI/blimpy/blob/master/CODE_OF_CONDUCT.md). Basically, be nice.
* If you have a small improvement, such as fixing a typo in the documentation, just open an issue and explain in as best detail as you can.
* If you have a moderate change that you've already coded up, issue a pull request with as much detail as possible about why it's awesome.
* We run continuous integration to check if new code breaks stuff, and checks code coverage. We probably won't merge code that breaks stuff, or lowers the coverage drastically. This probably means you'll need some unit tests. If you need help writing tests, please open a pull request anyway and we'll try to work with you to bring it up to a mergable state.
* If you issue a push request with new functionality, please write up some examples and unit tests so we can figure out what it does. 
* If you're planning major changes, please open up a discussion on the issue tracker first! We use `blimpy` in our telescope processing pipelines, so need to make sure our observation scripts won't break.
* `blimpy` is used in installations where Python 2.7 is still used. So code should, for the time being, be Py2 and Py3 compatible.

### Syle: PEP8 and docstrings

We favour the more concise [google docstrings](http://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings), 
and generally try to follow the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008/), but we do like to line up equal signs over multiple lines. 
That is: try and write PEP8 compliant code, but we don't mind well motivated deviations.

### Reporting issues or problems with the software 

* Please use the [github issues](https://github.com/UCBerkeleySETI/blimpy/issues) page to open an issue. Try and explain in the most detail, so we can recreate and fix it.

### Getting support

* [Open an issue](https://github.com/UCBerkeleySETI/blimpy/issues). If it's a bug report, or a feature request, please use the templates; otherwise try to be as descriptive as you can about what you're trying to do. Include as much relevant information as possible. 
