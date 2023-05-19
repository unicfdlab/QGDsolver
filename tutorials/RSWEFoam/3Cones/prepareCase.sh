funkySetFields -create -field b -expression 'max(max(max(0, 3.0 - 0.3*sqrt(sqr(pos().x - 10.0) + sqr(pos().y))), 1.0 - sqrt(sqr(pos().x + 7.5) + sqr(pos().y- 9.0))/8.0), 1.0 - sqrt(sqr(pos().x + 7.5) + sqr(pos().y+ 9.0))/8.0)' -dimension '[0 1 0 0 0 0 0]' -time 0

funkySetFields -create -field h -expression 'pos().x <= -21.5 ? 1.875 : 0.0' -dimension '[0 1 0 0 0 0 0]' -time 0