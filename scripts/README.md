NB: These scripts can generate massive temp files. On Windows, they're stored in your user profile by default. If you have a setup where your OS is on a little SSD and your data files are on a separate larger HDD, you may have trouble with the SSD filling up. In that case, change the default location of your temp directory with

    write("TMPDIR = '<your-desired-tempdir>'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

as per the advice at http://stackoverflow.com/questions/17107206/change-temporary-directory.
