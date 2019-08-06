%%quxianyaoshu(lo,R,a) %%quxianyaoshu(60,500,dms2degrees([28 36 20]))  
%%%%%%%%%%%加入桥墩偏距，直接搞桥墩偏距
disp('请输入半径（R），曲线长（lo)，起点坐标QD,终点坐标ED,交点坐标JD')
disp('格式如下：hhqxgx(R,lo,QD,ED,JD,lc)')
disp('比如：QD=[3273704.468,493811.52];ED=[3272419.04,495669.047];JD=[3273483.392,494887.186];R=607.2;lo=60;lc=80;')
disp('则输入以下命令行')
disp('hhqxgx(607.2,60,[3273704.468,493811.52],[3272419.04,495669.047],[3273483.392,494887.186],80')
%function hhqxgx(R,lo,QD,ED,JD,lc)
%输入里程；
close all
lc=80;
%请输入QD,JD和ED
QD=[3722.1,3830.5];
JD=[3501.7,4907];
ED=[2436.3,5688.9];
%请输入R,L0.
R=611.7;
lo=60;
%计算阿尔法角度,
[a1,s1]=cart2pol(JD(1)-QD(1),JD(2)-QD(2));
[a2,s2]=cart2pol(ED(1)-JD(1),ED(2)-JD(2));
a1=a1*180/pi;
a2=a2*180/pi;
a=-a1+a2;
%计算曲线要素
beitao=(lo/(2*R))*(180/pi);
detao=(1/3)*beitao;
m=0.5*lo-((lo^3)/(240*(R^2)));
p=(lo^2)/(24*R)-(lo^4)/(2688*(R^3));
T=m+(R+p)*tand(a/2);
L=2*lo+pi*R*(a-2*beitao)/180 ;
Eo=(R+p)*secd(a/2)-R;
q=2*T-L;
yo=lo^2/(6*R);
xo=lo-lo^3/(40*R^2);
%计算ZH点坐标
ZH=[QD(1)+(s1-T)*cosd(a1),QD(2)+(s1-T)*sind(a1)];
%计算HZ点坐标
HZ=[ED(1)-(s2-T)*cosd(a2),ED(2)-(s2-T)*sind(a2)];
%从ZH点计算到HY的坐标
kkkkk=1;

%计算起点里程
ys=rem(s1-T,1);%起点里程
for l=1-ys:1:60
%x=l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4);
%y=l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5);
XYzh2hy(kkkkk,1)=ZH(1)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*cosd(a1)-(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*sind(a1);
XYzh2hy(kkkkk,2)=ZH(2)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*sind(a1)+(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*cosd(a1);
QDzh2hy(kkkkk,1)=ZH(1)+((l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))+0.2*sind((1/R)*(l^2/(2*lo))))*cosd(a1)-(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5)-0.2*cosd((1/R)*(l^2/(2*lo))))*sind(a1);
QDzh2hy(kkkkk,2)=ZH(2)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4)+0.2*sind((1/R)*(l^2/(2*lo))))*sind(a1)+(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5)-0.2*cosd((1/R)*(l^2/(2*lo))))*cosd(a1);
kkkkk=kkkkk+1;
end
l=lo;
HY(1)=ZH(1)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*cosd(a1)-(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*sind(a1);
HY(2)=ZH(2)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*sind(a1)+(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*cosd(a1);
plot(XYzh2hy(:,1),XYzh2hy(:,2));
hold on;
plot(QDzh2hy(:,1),QDzh2hy(:,2),'k');



SR=(a-2*beitao)/180*pi*R;
%计算O点坐标
ak=a1-(180-a)/2;
OXY=JD-[(Eo+R)*cosd(ak),(Eo+R)*sind(ak)];
kkk=0;
%计算圆曲线中线
for aa=(ak-(a/2-beitao)+(1-ys)/R*180/pi):1/R*180/pi:(ak+(a/2-beitao))
kkk=kkk+1;
XYyqx(kkk,1)=OXY(1)+R*cosd(aa);
XYyqx(kkk,2)=OXY(2)+R*sind(aa);
%计算圆曲线上桥墩点位坐标
QDyqx(kkk,1)=OXY(1)+(R+0.4)*cosd(aa);
QDyqx(kkk,2)=OXY(2)+(R+0.4)*sind(aa);
end
plot(XYyqx(:,1),XYyqx(:,2));
plot(QDyqx(:,1),QDyqx(:,2),'k');

sssa=rem(SR-(1-ys),1);
%从HZ点计算到YH的坐标
KKKK=1;
for l=sssa:1:60
%x=l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4);
%y=l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5);

XYhz2yh(KKKK,1)=HZ(1)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*cosd(a2+180)+(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*sind(a2+180);
XYhz2yh(KKKK,2)=HZ(2)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*sind(a2+180)-(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*cosd(a2+180);
QDhz2yh(KKKK,1)=HZ(1)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4)+0.2*sind((1/R)*(l^2/(2*lo))))*cosd(a2+180)+(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5)-0.2*cosd((1/R)*(l^2/(2*lo))))*sind(a2+180);
QDhz2yh(KKKK,2)=HZ(2)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4)+0.2*sind((1/R)*(l^2/(2*lo))))*sind(a2+180)-(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5)-0.2*cosd((1/R)*(l^2/(2*lo))))*cosd(a2+180);
KKKK=KKKK+1;
end
l=lo;
YH(1)=HZ(1)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*cosd(a2+180)+(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*sind(a2+180);
YH(2)=HZ(2)+(l-l^5/(40*R^2*lo^2)+l^9/(3456*R^4*lo^4))*sind(a2+180)-(l^3/(6*R*lo)-l^7/(336*R^3*lo^3)+l^11/(42240*R^5*lo^5))*cosd(a2+180);
plot(XYhz2yh(:,1),XYhz2yh(:,2));
plot(QDhz2yh(:,1),QDhz2yh(:,2),'k');
%从HY点计算到YH的坐标








      disp('注意事项如下，请在代码中修改你的R（每个人，每个小组都不一样），lo(一样），里程（一样）');
      disp('--------------------------------------------------------------------------------------------------------');
    disp('黑色为墩台线,QD**2**,表示**到**这段中线的桥墩坐标');
    disp('黑色为墩台线,QD**2**,表示**到**这段中线的桥墩坐标');
    disp('黑色为墩台线,QD**2**,表示**到**这段中线的桥墩坐标');
    disp('--------------------------------------------------------------------------------------------------------');
    disp(' XYzh2hy为ZH到HY的坐标点');
    disp(' XYhz2yh为HZ到YH的坐标点');
    disp(' XYyqx为HY到YH的坐标点');
        disp('鑫哥');
    text(QD(1),QD(2),'QD','horiz','left','color','r','fontsize',10);
     text(JD(1),JD(2),'JD','horiz','left','color','r','fontsize',10);
   text(ED(1),ED(2),'ED','horiz','left','color','r','fontsize',10);
     text(ZH(1),ZH(2),'ZH','horiz','left','color','r','fontsize',10);
     text(HZ(1),HZ(2),'HZ','horiz','left','color','r','fontsize',10);
         text(YH(1),YH(2),'YH','horiz','left','color','r','fontsize',10);
     text(HY(1),HY(2),'HY','horiz','left','color','r','fontsize',10);

     %载入导线坐标图;
     dx=[3273442.18000000000000000 	494585.17700000000000000 
	
3273448.87724612000000000 	494656.49984010800000000 
	
3273444.21890341000000000 	494723.40269034100000000 
	
3273436.90359182000000000 	494788.11444492600000000 
	
3273392.59355999000000000 	494852.05766109300000000 
	
3273348.43925840000000000 	494912.28817138200000000 
	
3273349.54877251000000000 	494985.45635436700000000 
	
3273388.06168143000000000 	495045.38349452000000000 
	
3273387.27473207000000000 	495097.59376429600000000 
	
3273336.99901743000000000 	495100.97202560200000000 
	
3273284.92203112000000000 	495097.20150399400000000 
	
3273258.06514289000000000 	495154.24691697800000000 
	
3273184.91549385000000000 	495150.66031005400000000 
	
3273115.34463699000000000 	495161.17597306700000000 
	
3273024.48723111000000000 	495156.73787290400000000 
	
3272950.97948669000000000 	495110.11443847100000000 
	
3272884.18373020000000000 	495079.05409712900000000 
	
3272838.68200001000000000 	495110.61900000100000000 
	
];
    %----------------------------------------------------------------
%plot(dx(:,1),dx(:,2),'r');
%plot(dx(:,1),dx(:,2),'o');
  %text(dx(10,1),dx(10,2)+15,'9','horiz','left','color','r','fontsize',10);
  %text(dx(9,1),dx(9,2)+15,'8','horiz','left','color','r','fontsize',10);
 %text(dx(11,1),dx(11,2)+15,'10','horiz','left','color','r','fontsize',10);
   %text(dx(7,1),dx(7,2)+15,'6','horiz','left','color','r','fontsize',10);
 %text(dx(8,1),dx(8,2)+15,'7','horiz','left','color','r','fontsize',10);
     %text(3273361.9700,495034.8210+15,'JM07','horiz','left','color','r','fontsize',10);%加密点
   % plot(3273361.9700,495034.8210+15,'o');
    %JM=[3273361.9700,495034.8210];
    
    %----------------------------------------------------------------
        %删除上一次的excel文件
        delete('坐标集合.xls');
     %导出EXCEL代码
     xlswrite('坐标集合', QDzh2hy, '墩点位 ZH到HY,里程0-60');
    xlswrite('坐标集合', QDyqx, ' 桥墩点位 HY到YH,里程60-R乘以圆心角');
     xlswrite('坐标集合', QDhz2yh, '桥墩点位 YH到HZ,里程');
   xlswrite('坐标集合', XYzh2hy, '中线坐标 ZH到HY,里程0-60');
    xlswrite('坐标集合', XYyqx, '中线坐标 HY到YH,里程60-385');
     xlswrite('坐标集合', XYhz2yh, '中线坐标 HZ到YH,里程');
     disp('第一页 桥墩点位 ZH到HY;第二页桥墩点位 HY到YH,第三页桥墩点位 YH到HZ,');
          disp('第四页 中线坐标 ZH到HY;第四页中线坐标 HY到YH,第五页中线坐标 YH到HZ,');
          %偏距检测
          for i=1:60
              dhzyh(i)=norm(QDhz2yh(i,:)-XYhz2yh(i,:));
          end
          for i=1:60
              dzhhy(i)=norm(QDzh2hy(i,:)-XYzh2hy(i,:));
          end
          for i=1:386
              dyqx(i)=norm(QDyqx(i,:)-XYyqx(i,:));
          end
          dhzyh=dhzyh';
          dzhhy=dzhhy';
          dyqx=dyqx';
          
          
          %加入文字说明
          title('导线与曲线位置示意图');
          
          %里程检测：
          norm(XYhz2yh(60,:)-XYyqx(386,:))
          norm(XYzh2hy(60,:)-XYyqx(1,:))
          %norm(XYzh2hy(1,:)-ZH)
          norm(XYzh2hy(1,:)-QD)-(s1-T)+ys
          for i=1:59
              lchzyh(i)=norm(XYhz2yh(i+1,:)-XYhz2yh(i,:));
          end
          for i=1:59
              lchhy(i)=norm(XYzh2hy(i+1,:)-XYzh2hy(i,:));
          end
          for i=1:385
              lcyqx(i)=norm(XYyqx(i+1,:)-XYyqx(i,:));
          end
          disp('ZH点坐标里程如下')
          norm(XYzh2hy(1,:)-QD)+lc
