if hash libtoolize 2>&-
then
    libtoolize && aclocal && automake --add-missing && autoconf
else
    glibtoolize && aclocal && automake --add-missing && autoconf
fi
