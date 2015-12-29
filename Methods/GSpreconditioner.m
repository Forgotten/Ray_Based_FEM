function uPrecond = GSpreconditioner(MArray, indx0, indx1, indxn ,indxnp, ...
                                     fIntLocal )



localLim = [0 cumsum(localSizes)];

%% local solve + downwards sweep
for ii = 1:nSub
   
        % making a local copy of the local rhs
    rhsLocaltemp = fIntLocal{ii};

    if ii ~=1
        rhsLocaltemp(indx0{ii}) =  rhsLocaltemp(indx0{ii}) +...
                   MArray{ii}.H(indx0{ii},indx1{ii})*u_np((ii-2)*n + index);
        rhsLocaltemp(indx1{ii}) = rhsLocaltemp(indx1{ii}) -...
                    MArray{ii}.H(indx1{ii},indx0{ii})*u_n((ii-2)*n + index);
    end

        % solving the rhs
        vDown = MArray{ii}.solveInt(rhsLocaltemp);

        % extracting the traces
        if ii ~= nSub
            u_n((ii-1)*n  + index) = vDown(indxn{ii});
            u_np((ii-1)*n + index) = vDown(indxnp{ii});
        end
end

%% upwards sweep + reflections + reconstruction

uPrecond = zeros(size(fInt,1),1);

     for ii = nSub:-1:1
      
        % making a copy of the parititioned source
         rhsLocaltemp = fIntLocal{ii};

        % adding the source at the boundaries
        if ii~= 1
            % we need to be carefull at the edges
            rhsLocaltemp(indx1{ii})  = rhsLocaltemp(indx1{ii}) - ...
               MArray{ii}.H(indx1{ii},indx0{ii})*u_n((ii-2)*n + index);
            rhsLocaltemp(indx0{ii})  =  rhsLocaltemp(indx0{ii}) ...
                + MArray{ii}.H(indx0{ii},indx1{ii})*u_np((ii-2)*n + index);
        end
        if ii~= nSub
            rhsLocaltemp(indxnp{ii}) = rhsLocaltemp(indxnp{ii}) + ...
                MArray{ii}.H(indxnp{ii},indxn{ii})*u_0((ii)*n + index);
            rhsLocaltemp(indxn{ii})  = rhsLocaltemp(indxn{ii}) - ...
                MArray{ii}.H(indxn{ii},indxnp{ii})*u_1((ii)*n + index);
        end
 
        % solving the local problem
        uUp = MArray{ii}.solveInt(rhsLocaltemp);

        if ii > 1
            u_0((ii-1)*n + index) = uUp(indx0{ii});
            u_1((ii-1)*n + index) = uUp(indx1{ii}) - u_np((ii-2)*n + index);
        end

        % reconstructing the problem on the fly
        uPrecond(localLim(ii)+1:localLim(ii+1)) = uUp(indIntLocal{ii});
    end