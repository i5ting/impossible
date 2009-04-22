PROG=		bi
CFLAGS=		-g -I/usr/local/include
LDFLAGS=	-lm -L/usr/local/lib -lsndfile

.include <bsd.prog.mk>
