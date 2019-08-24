### Test Data

------
To remove unwanted bulk from this repository,
the test data will not be kept in the repo, but will be stored instead on a public server.

Currently, they are [here](http://blpd0.ssl.berkeley.edu/Voyager_data/).

If more data needs to be added, please modify `download_data.sh`.

e.g.
```bash
curl --url "http://blpd0.ssl.berkeley.edu/Voyager_data/file_name" -o ./test_data/file_name
```