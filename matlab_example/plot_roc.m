% predict-�������Բ��Լ��ķ�����
% ground. _truth -���Լ�����ȷ��ǩ,����ֻ���Ƕ�����,��0��1
% auc-����ROC���ߵ����� �µ����
function auc = plot_roc( predict, ground_truth)
%��ʼ��Ϊ( 1.0, 1.0 )
%�����ground_ truth������������Ŀpos_ num�͸���������Ŀneg. num
pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==0);
m=size(ground_truth,1);
[pre,Index]= sort(predict);
ground_truth=ground_truth(Index);
x=zeros(m+1,1);
y=zeros(m+1,1);
auc=0;
x(1)=1;
y(1)=1; 
for i=2:m
TP=sum(ground_truth(i:m) == 1);
FP= sum(ground_truth(i:m) == 0);
x(i)=FP/neg_num;
y(i)=TP/pos_num;
auc=auc + (y(i)+y(i-1))*(x(i-1)-x(i))/2;
end;
x(m+ 1)=0;y(m+1)=0;
auc= auc+y(m)*x(m)/2;
% plot(x,y);

plot(x,y,'-bo','LineWidth',2,'MarkerSize',3);
% xlabel('�鱨����');
% ylabel('���и���');
% title('ROC����ͼ');

end

%%
% function  auc = plot_roc( predict, ground_truth )
% % INPUTS
% %  predict       - �������Բ��Լ��ķ�����
% %  ground_truth - ���Լ�����ȷ��ǩ,����ֻ���Ƕ����࣬��0��1
% % OUTPUTS
% %  auc            - ����ROC���ߵ������µ����
% 
% %��ʼ��Ϊ��1.0, 1.0��
% x = 1.0;
% y = 1.0;
% %�����ground_truth������������Ŀpos_num�͸���������Ŀneg_num
% pos_num = sum(ground_truth==1);
% neg_num = sum(ground_truth==0);
% %���ݸ���Ŀ���Լ������x�����y��Ĳ���
% x_step = 1.0/neg_num;
% y_step = 1.0/pos_num;
% %���ȶ�predict�еķ��������ֵ���մ�С��������
% [predict,index] = sort(predict);
% ground_truth = ground_truth(index);
% %��predict�е�ÿ�������ֱ��ж�������FP������TP
% %����ground_truth��Ԫ�أ�
% %��ground_truth[i]=1,��TP������1����y�᷽���½�y_step
% %��ground_truth[i]=0,��FP������1����x�᷽���½�x_step
% for i=1:length(ground_truth)
%     if ground_truth(i) == 1
%         y = y - y_step;
%     else
%         x = x - x_step;
%     end
%     X(i)=x;
%     Y(i)=y;
% end
% %����ͼ��     
% plot(X,Y,'-ro','LineWidth',2,'MarkerSize',3);
% xlabel('�鱨����');
% ylabel('���и���');
% title('ROC����ͼ');
% %����С���ε����,����auc
% auc = -trapz(X,Y);  
% 
% end

% https://blog.csdn.net/xmu_jupiter/article/details/21885299