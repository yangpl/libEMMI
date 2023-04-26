cp rho11_init rho11_init2
printf "\x0" | dd of=rho11_init2 bs=1 seek=0 count=1 conv=notrunc
printf "\x1" | dd of=rho11_init2 bs=1 seek=0 skip=8 count=1 conv=notrunc
