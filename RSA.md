## RSA Starter 1
### Solution
Python have `pow()` function

### Code

```python
print(pow(101,17,22663))
```

## RSA Starter 2
### Solution
Print $12^{p \times q} \mod e$

### Code

```python
e = 65537
p = 17
q = 23
n = 17*23
print(pow(12, e, n))
```

## RSA Starter 3
### Solution
Because $N = p \times q$ and $p$, $q$ are primes so
$$
\phi\left({n}\right) = \left({p - 1}\right) \times \left({q - 1}\right)
$$

### Code

```python
print(857504083339712752489993810776 * 1029224947942998075080348647218)
```

## RSA Starter 4
### Solution
$$
d = e^{-1} \mod \phi\left({n}\right)
$$

### Code

```python
from Crypto.Util.number import inverse
e = 65537
p = 857504083339712752489993810777
q = 1029224947942998075080348647219
phi = (p - 1) * (q - 1)
d = inverse(e, phi)
print(d)
```

## RSA Starter 5
### Solution
First, I used [factordb](http://factordb.com/) to find the factors of $n$
![[RSAStarter5.png]]
We have $p$ and $q$ now.
Calculate in sequence:
$$
\begin{align*}
\phi  &= \left( {p - 1} \right) \times \left( {q - 1} \right)\\
d &= {e^{ - 1}}\bmod \phi \\
m &= {c^d}\bmod n
\end{align*}
$$

### Code

```python
from Crypto.Util.number import inverse
e = 65537
p = 857504083339712752489993810777
q = 1029224947942998075080348647219
n = p * q
phi = (p - 1) * (q - 1)
d = inverse(e, phi)
c = 77578995801157823671636298847186723593814843845525223303932
print(pow(c, d, n))
```

## RSA Starter 6
### Solution
SHA256 first, and then encrypt
$$
c = m^{d} \mod n
$$

### Code

```python
from Crypto.Util.number import inverse, bytes_to_long
from hashlib import sha256
N = 15216583654836731327639981224133918855895948374072384050848479908982286890731769486609085918857664046075375253168955058743185664390273058074450390236774324903305663479046566232967297765731625328029814055635316002591227570271271445226094919864475407884459980489638001092788574811554149774028950310695112688723853763743238753349782508121985338746755237819373178699343135091783992299561827389745132880022259873387524273298850340648779897909381979714026837172003953221052431217940632552930880000919436507245150726543040714721553361063311954285289857582079880295199632757829525723874753306371990452491305564061051059885803
d = 11175901210643014262548222473449533091378848269490518850474399681690547281665059317155831692300453197335735728459259392366823302405685389586883670043744683993709123180805154631088513521456979317628012721881537154107239389466063136007337120599915456659758559300673444689263854921332185562706707573660658164991098457874495054854491474065039621922972671588299315846306069845169959451250821044417886630346229021305410340100401530146135418806544340908355106582089082980533651095594192031411679866134256418292249592135441145384466261279428795408721990564658703903787956958168449841491667690491585550160457893350536334242689

s = "crypto{Immut4ble_m3ssag1ng}"
sh = sha256(s.encode()).digest()
m = pow(bytes_to_long(sh), d, N)
print(m)
```

## Factoring
### Solution
Yes, factordb
![[Factoring.png]]

## Inferius Prime
### Solution
Use factordb again and you will get $p$ and $q$
![[Inferiusprime.png]]

After that, decrypt like RSA Starter 5

### Code

```python
from Crypto.Util.number import inverse, bytes_to_long, long_to_bytes

n = 742449129124467073921545687640895127535705902454369756401331
e = 3
ct = 39207274348578481322317340648475596807303160111338236677373

p = 752708788837165590355094155871
q = 986369682585281993933185289261

phi = (p -1) * (q - 1)
d = inverse(e, phi)
m = pow(ct, d, n)
print(long_to_bytes(m))
```

## Monoprime
### Solution
If $n$ is prime, $\phi\left({n}\right) = n - 1$

### Code

```python
from Crypto.Util.number import inverse, bytes_to_long, long_to_bytes
n = 171731371218065444125482536302245915415603318380280392385291836472299752747934607246477508507827284075763910264995326010251268493630501989810855418416643352631102434317900028697993224868629935657273062472544675693365930943308086634291936846505861203914449338007760990051788980485462592823446469606824421932591                                                                  
e = 65537
ct = 161367550346730604451454756189028938964941280347662098798775466019463375610700074840105776873791605070092554650190486030367121011578171525759600774739890458414593857709994072516290998135846956596662071379067305011746842247628316996977338024343628757374524136260758515864509435302781735938531030576289086798942  
phi = n - 1
d = inverse(e, phi)
m = pow(ct, d, n)
print(long_to_bytes(m))
```

## Square Eyes
### Solution
You can use factordb to factorize first and we have this formula:
$$
\phi \left( n \right) = n\left( {1 - \frac{1}{p}} \right) = n - p = p\left( {p - 1} \right)
$$

### Code

```python
from Crypto.Util.number import inverse, bytes_to_long, long_to_bytes
import math

n = 535860808044009550029177135708168016201451343147313565371014459027743491739422885443084705720731409713775527993719682583669164873806842043288439828071789970694759080842162253955259590552283047728782812946845160334801782088068154453021936721710269050985805054692096738777321796153384024897615594493453068138341203673749514094546000253631902991617197847584519694152122765406982133526594928685232381934742152195861380221224370858128736975959176861651044370378539093990198336298572944512738570839396588590096813217791191895941380464803377602779240663133834952329316862399581950590588006371221334128215409197603236942597674756728212232134056562716399155080108881105952768189193728827484667349378091100068224404684701674782399200373192433062767622841264055426035349769018117299620554803902490432339600566432246795818167460916180647394169157647245603555692735630862148715428791242764799469896924753470539857080767170052783918273180304835318388177089674231640910337743789750979216202573226794240332797892868276309400253925932223895530714169648116569013581643192341931800785254715083294526325980247219218364118877864892068185905587410977152737936310734712276956663192182487672474651103240004173381041237906849437490609652395748868434296753449
e = 65537
ct = 222502885974182429500948389840563415291534726891354573907329512556439632810921927905220486727807436668035929302442754225952786602492250448020341217733646472982286222338860566076161977786095675944552232391481278782019346283900959677167026636830252067048759720251671811058647569724495547940966885025629807079171218371644528053562232396674283745310132242492367274184667845174514466834132589971388067076980563188513333661165819462428837210575342101036356974189393390097403614434491507672459254969638032776897417674577487775755539964915035731988499983726435005007850876000232292458554577437739427313453671492956668188219600633325930981748162455965093222648173134777571527681591366164711307355510889316052064146089646772869610726671696699221157985834325663661400034831442431209123478778078255846830522226390964119818784903330200488705212765569163495571851459355520398928214206285080883954881888668509262455490889283862560453598662919522224935145694435885396500780651530829377030371611921181207362217397805303962112100190783763061909945889717878397740711340114311597934724670601992737526668932871436226135393872881664511222789565256059138002651403875484920711316522536260604255269532161594824301047729082877262812899724246757871448545439896
p = 23148667521998097720857168827790771337662483716348435477360567409355026169165934446949809664595523770853897203103759106983985113264049057416908191166720008503275951625738975666019029172377653170602440373579593292576530667773951407647222757756437867216095193174201323278896027294517792607881861855264600525772460745259440301156930943255240915685718552334192230264780355799179037816026330705422484000086542362084006958158550346395941862383925942033730030004606360308379776255436206440529441711859246811586652746028418496020145441513037535475380962562108920699929022900677901988508936509354385660735694568216631382653107
phi = (p - 1) * (p)
d = inverse(e, phi)
m = pow(ct, d, n)
print(long_to_bytes(m))
```

## Manyprime
### Solution
Factordb again, it's too strong
![[Manyprime.png]]

### Code

```python
from Crypto.Util.number import inverse, bytes_to_long, long_to_bytes
import requests


n = 580642391898843192929563856870897799650883152718761762932292482252152591279871421569162037190419036435041797739880389529593674485555792234900969402019055601781662044515999210032698275981631376651117318677368742867687180140048715627160641771118040372573575479330830092989800730105573700557717146251860588802509310534792310748898504394966263819959963273509119791037525504422606634640173277598774814099540555569257179715908642917355365791447508751401889724095964924513196281345665480688029639999472649549163147599540142367575413885729653166517595719991872223011969856259344396899748662101941230745601719730556631637
e = 65537
ct = 320721490534624434149993723527322977960556510750628354856260732098109692581338409999983376131354918370047625150454728718467998870322344980985635149656977787964380651868131740312053755501594999166365821315043312308622388016666802478485476059625888033017198083472976011719998333985531756978678758897472845358167730221506573817798467100023754709109274265835201757369829744113233607359526441007577850111228850004361838028842815813724076511058179239339760639518034583306154826603816927757236549096339501503316601078891287408682099750164720032975016814187899399273719181407940397071512493967454225665490162619270814464

arr = [9282105380008121879, 9303850685953812323, 9389357739583927789, 10336650220878499841, 10638241655447339831, 11282698189561966721, 11328768673634243077, 11403460639036243901, 11473665579512371723, 11492065299277279799, 11530534813954192171, 11665347949879312361, 12132158321859677597, 12834461276877415051, 12955403765595949597, 12973972336777979701, 13099895578757581201, 13572286589428162097, 14100640260554622013, 14178869592193599187, 14278240802299816541, 14523070016044624039, 14963354250199553339, 15364597561881860737, 15669758663523555763, 15824122791679574573, 15998365463074268941, 16656402470578844539, 16898740504023346457, 17138336856793050757, 17174065872156629921, 17281246625998849649] 

phi = 1
for i in arr:
    phi *= (i - 1)
d = inverse(e, phi)
m = pow(ct, d, n)
print(long_to_bytes(m))
```

## Salty
### Solution
I figured it out that there is a [library](https://github.com/ryosan-470/factordb-python) to use factordb in python.
Just use it and decrypt.

### Code

```python
from Crypto.Util.number import inverse, bytes_to_long, long_to_bytes
from factordb.factordb import FactorDB

n = 110581795715958566206600392161360212579669637391437097703685154237017351570464767725324182051199901920318211290404777259728923614917211291562555864753005179326101890427669819834642007924406862482343614488768256951616086287044725034412802176312273081322195866046098595306261781788276570920467840172004530873767                                                                  
e = 1
ct = 44981230718212183604274785925793145442655465025264554046028251311164494127485

f = FactorDB(n)
f.connect()

arr = f.get_factor_list()

phi = 1
for i in arr:
    phi *= (i - 1)
d = inverse(e, phi)
m = pow(ct, d, n)
print(long_to_bytes(m))

```

## 