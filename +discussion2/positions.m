function [position,velocity,horizon] = positions(cfg,bank,kind)
%POSITIONS Original ROE entries or uniform-volume stationary probes.
arguments
    cfg (1,1) struct
    bank (1,1) struct
    kind (1,1) string {mustBeMember(kind,["dynamic","static"])}
end
u = bank.u;
azimuth = 2*pi*u(:,1);
cosPhi = u(:,2);
direction = [sqrt(1-cosPhi.^2).*cos(azimuth), ...
    sqrt(1-cosPhi.^2).*sin(azimuth),cosPhi];
position = cfg.system.maxRange_m*direction;
if kind == "static"
    position = position.*nthroot(max(u(:,8),realmin),3);
    velocity = zeros(bank.N,3);
    horizon = inf(bank.N,1);
    return
end
velocityAz = 2*pi*u(:,3);
velocityCos = 2*u(:,4)-1;
unitVelocity = [sqrt(1-velocityCos.^2).*cos(velocityAz), ...
    sqrt(1-velocityCos.^2).*sin(velocityAz),velocityCos];
outward = sum(unitVelocity.*direction,2)>0;
unitVelocity(outward,:) = -unitVelocity(outward,:);
velocity = cfg.system.targetSpeed_m_s*unitVelocity;
sphereTime = -2*sum(position.*velocity,2)./sum(velocity.^2,2);
groundTime = inf(bank.N,1);
down = velocity(:,3)<0;
groundTime(down) = -position(down,3)./velocity(down,3);
horizon = max(0,min([sphereTime,groundTime, ...
    repmat(cfg.maxTime_s,bank.N,1)],[],2));
end
