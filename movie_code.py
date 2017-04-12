    plot(dump[:, 0])
    ylim([0.9899, 1.0101])
    xlim([0, 256])
    title('$\\rho$')
    savefig('rho' + '%04d'%i + '.png')
    clf()

ffmpeg -f image2 -i rho%04d.png -vcodec mpeg4 -mbd rd -trellis 2 -cmp 2 -g 300 -pass 1 -r 1 -b 18000000 movie.mp4
